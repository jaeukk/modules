"""Run bounded, serial before/after contact-reduction experiments with preserved receipts."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import numpy as np

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
VARIANTS={'one_2t':(500,1,2,299),'one_4t':(500,1,4,299),'five_4t':(2000,5,4,0)}
SEEDS=[20260928,20260929,20260930]


def digest(path):
    """Return the SHA-256 of an evidence file."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path,value):
    """Atomically checkpoint JSON evidence."""
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)


def environment(threads):
    """Match both the requested OpenMP team and runtime limit to physical CPUs."""
    cpus=[0,2,4,6][:threads]
    os.sched_setaffinity(0,cpus)
    return dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_THREAD_LIMIT=str(threads),
                OMP_DYNAMIC='false',OMP_WAIT_POLICY='PASSIVE',OMP_PROC_BIND='false',
                OPENBLAS_NUM_THREADS=str(threads),MKL_NUM_THREADS=str(threads))


def compare(a,b,n):
    """Compare fixed particle IDs with periodic wrapping and optimal translation."""
    x=np.loadtxt(a,skiprows=1);y=np.loadtxt(b,skiprows=1)
    assert len(x)==len(y)==n and np.array_equal(x[:,[0,1,4]],y[:,[0,1,4]])
    delta=y[:,2:4]-x[:,2:4];delta-=np.floor(delta+.5)
    for d in range(2):
        vals=np.sort(delta[:,d]);best=(float('inf'),0.)
        for cut in range(n):
            mean=(vals.sum()+cut)/n
            residual=vals-mean;residual-=np.floor(residual+.5)
            score=float(residual@residual)
            if score<best[0]:best=(score,mean)
        delta[:,d]-=best[1]
    delta-=np.floor(delta+.5)
    distance=np.linalg.norm(delta,axis=1)*np.sqrt(n)
    return dict(bitwise_identical=Path(a).read_bytes()==Path(b).read_bytes(),
                mean=float(distance.mean()),rms=float(np.sqrt(np.mean(distance**2))),
                maximum=float(distance.max()),mean_relative=float(distance.mean()/np.sqrt(n)))


def report(root):
    """Summarize available phases without conflating repeat and backend differences."""
    lines=['# Contact-repulsion and fixed-order summation study','',
           'Three matched seeds, two independent process repeats per backend. All runs use '
           'phi_total=0.2, equal radii sqrt(0.2/pi), density=1, u0=e0=1, and the unchanged '
           'three-pass L-BFGS routine (2000 evaluations/pass, 600 seconds/configuration). '
           'One-component N500 uses the previous 299 modes (chi=299/998); five-component '
           'N2000 uses 80 complete-shell modes per species (chi=80/798).', '',
           'Times are optimization wall times. Each endpoint is independently checked by direct '
           'Fourier sums and all-pair contact forces. Displacements use fixed labels, periodic '
           'wrapping, and optimal translation. Acceptance is E<=1e-11, maxS<=1e-10, '
           'max force component<=1e-7, max relative overlap<=1e-6.','']
    comparisons=[]
    for phase in ['baseline','fixed']:
        phase_dir=root/phase
        if not phase_dir.exists():continue
        lines += [f'## {phase}', '', '| Case | Seed | Vanilla seconds (r0/r1) | FINUFFT seconds (r0/r1) | Vanilla repeat mean | FINUFFT repeat mean | Cross-backend mean (r0/r1) | All accepted |',
                  '|---|---|---|---|---:|---:|---|---|']
        for variant,(n,species,threads,modes) in VARIANTS.items():
            for seed in SEEDS:
                paths={(b,r):phase_dir/f'{variant}_{seed}_{b}_{r}' for b in ['direct','finufft'] for r in range(2)}
                if not all(p.with_suffix('.json').exists() for p in paths.values()):continue
                rows={k:json.loads(p.with_suffix('.json').read_text()) for k,p in paths.items()}
                assert len({digest(p/'initial.tsv') for p in paths.values()})==1
                assert len({digest(p/'modes.tsv') for p in paths.values()})==1
                same={b:compare(paths[b,0]/'final.tsv',paths[b,1]/'final.tsv',n) for b in ['direct','finufft']}
                cross=[compare(paths['direct',r]/'final.tsv',paths['finufft',r]/'final.tsv',n) for r in range(2)]
                comparisons.append(dict(phase=phase,variant=variant,seed=seed,same_backend=same,cross_backend=cross))
                times={b:'/'.join(f"{rows[b,r]['optimization_seconds']:.3f}" for r in range(2)) for b in ['direct','finufft']}
                accepted=all(v['direct_reference']['ground_state'] for v in rows.values())
                lines.append(f"| {variant} | {seed} | {times['direct']} | {times['finufft']} | {same['direct']['mean']:.6e} | {same['finufft']['mean']:.6e} | {cross[0]['mean']:.6e}/{cross[1]['mean']:.6e} | {accepted} |")
    save(root/'comparisons.json',comparisons)
    lines += ['', 'Per-run JSON files retain energies, forces, maximum overlap, optimizer statuses and '
              'evaluation counts, transform counters, timing, and input hashes. Probe JSON files retain '
              '100 repeated fixed-coordinate evaluations plus independent-force and finite-difference checks. '
              'The baseline and fixed executables and their build commands/source hashes are preserved separately.', '']
    (root/'REPORT.md').write_text('\n'.join(lines))


def main():
    """Build against unchanged core objects, then execute only this finite test matrix."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase',choices=['baseline','fixed'])
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--reference-build',type=Path,required=True)
    args=parser.parse_args();root=args.output.resolve();root.mkdir(exist_ok=True,parents=True)
    phase=root/args.phase;phase.mkdir(exist_ok=False)
    original=args.reference_build.resolve()
    commands=json.loads((original/'commands.json').read_text())
    original_hashes=json.loads((original/'source_hashes.json').read_text())
    for source,sha in original_hashes.items():
        if Path(source).name not in ['experiment.cpp','MultiStealthyPotential.cpp']:
            assert digest(source)==sha, f'Changed shared source: {source}'
    common=commands[0][:commands[0].index('-c')]
    contact=REPO/'apps/CCOptimization/MultiStealthyPotential.cpp'
    contact_obj=phase/'MultiStealthyPotential.o'
    compile_contact=common+['-c',str(contact),'-o',str(contact_obj)]
    subprocess.run(compile_contact,check=True,capture_output=True)
    save(phase/'build.json',dict(contact_sha256=digest(contact),contact_command=compile_contact,
         harness_sha256=digest(HERE/'experiment.cpp'),shared_harness_sha256=digest(HERE.parent/'five_component/experiment.cpp'),
         reused_object_hashes={x:digest(x) for x in commands[-1] if x.endswith('.o') and not x.endswith(('/experiment.o','/MultiStealthyPotential.o'))}))
    for variant,(n,species,threads,modes) in VARIANTS.items():
        env=environment(threads);binary=phase/variant;obj=phase/f'{variant}.o'
        compile_cmd=common+[f'-DCC_TEST_N={n}',f'-DCC_TEST_SPECIES={species}',f'-DCC_TEST_THREADS={threads}',f'-DCC_TEST_MODES={modes}',
                            '-c',str(HERE/'experiment.cpp'),'-o',str(obj)]
        subprocess.run(compile_cmd,check=True,capture_output=True)
        link=[str(obj) if x.endswith('/experiment.o') else str(contact_obj) if x.endswith('/MultiStealthyPotential.o') else x for x in commands[-1]]
        link[link.index('-o')+1]=str(binary);subprocess.run(link,check=True,capture_output=True)
        save(phase/f'{variant}_build.json',dict(commands=[compile_cmd,link],binary_sha256=digest(binary),
             thread_environment={k:v for k,v in env.items() if k.startswith(('OMP_','OPENBLAS_','MKL_'))},cpu_affinity=sorted(os.sched_getaffinity(0))))
        probe=subprocess.run([str(binary),'probe'],env=env,capture_output=True,text=True,check=True,timeout=120)
        save(phase/f'{variant}_probe.json',json.loads(probe.stdout));print(args.phase,variant,'probe passed',flush=True)
        for seed in SEEDS:
            for repeat in range(2):
                for backend in (['direct','finufft'] if (seed+repeat)%2==0 else ['finufft','direct']):
                    dest=phase/f'{variant}_{seed}_{backend}_{repeat}'
                    cmd=[str(binary),str(dest),str(seed),'0.3' if species==1 else '0.1',backend]
                    with dest.with_suffix('.stdout').open('w') as out,dest.with_suffix('.stderr').open('w') as err:
                        subprocess.run(cmd,env=env,stdout=out,stderr=err,check=True,timeout=660)
                    row=json.loads(dest.with_suffix('.stdout').read_text())
                    row.update(command=cmd,initial_sha256=digest(dest/'initial.tsv'),modes_sha256=digest(dest/'modes.tsv'))
                    save(dest.with_suffix('.json'),row)
                    print(args.phase,dest.name,round(row['optimization_seconds'],3),row['direct_reference']['ground_state'],flush=True)
            report(root)
    report(root)


if __name__=='__main__':main()
