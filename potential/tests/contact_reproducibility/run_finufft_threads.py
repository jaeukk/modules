"""Isolate FINUFFT nthreads=1 using the preserved four-thread fixed-sum experiment."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess

from run import compare, save

HERE=Path(__file__).resolve().parent


def digest(path):
    """Hash an input, object, or executable for the audit receipt."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit_log(text):
    """Require proof that actual FINUFFT calls used the intended override."""
    match=re.search(r'FINUFFT_THREAD_CONTROL nthreads=(\d+) calls=(\d+) unexpected_incoming=(\d+) nonzero_status=(\d+)',text)
    if not match:raise RuntimeError('missing FINUFFT thread-control receipt')
    result=dict(zip(['nthreads','calls','unexpected_incoming','nonzero_status'],map(int,match.groups())))
    assert result['nthreads']==1 and result['unexpected_incoming']==result['nonzero_status']==0,result
    return result


def main():
    """Execute twelve isolated matched-state optimizations without changing production settings."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixed-results',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();old=args.fixed_results.resolve();root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    build=json.loads((old/'five_4t_build.json').read_text());commands=build['commands']
    assert digest(old/'five_4t')==build['binary_sha256']
    source=HERE/'finufft_threads.cpp';obj=root/'finufft_threads.o';binary=root/'five_omp4_finufft1'
    common=commands[0][:commands[0].index('-c')]
    compile_cmd=common+['-c',str(source),'-o',str(obj)]
    subprocess.run(compile_cmd,check=True,capture_output=True)
    link=commands[1].copy();link[link.index('-o')+1]=str(binary)
    link.extend([str(obj),'-Wl,--wrap=finufft2d1,--wrap=finufft_makeplan,--wrap=finufft_setpts,--wrap=finufft_execute'])
    subprocess.run(link,check=True,capture_output=True)
    os.sched_setaffinity(0,build['cpu_affinity'])
    env=dict(os.environ,**build['thread_environment'])
    assert env['OMP_NUM_THREADS']==env['OMP_THREAD_LIMIT']=='4'
    save(root/'build.json',dict(commands=[compile_cmd,link],source_sha256=digest(source),
         reused_objects={p:digest(p) for p in commands[1] if p.endswith('.o')},
         binary_sha256=digest(binary),cpu_affinity=build['cpu_affinity'],thread_environment=build['thread_environment']))
    proc=subprocess.run([str(binary),'probe'],env=env,capture_output=True,text=True,check=True,timeout=120)
    (root/'probe.stderr').write_text(proc.stderr)
    probe=json.loads(proc.stdout);probe['thread_control']=audit_log(proc.stderr)
    save(root/'probe.json',probe)
    print('Probe:',json.dumps(probe),flush=True)
    rows=[];pairs=[]
    for seed in [20260928,20260929,20260930]:
        for repeat in range(2):
            for backend in (['direct','finufft'] if (seed+repeat)%2==0 else ['finufft','direct']):
                name=f'{seed}_{backend}_{repeat}';dest=root/name
                cmd=[str(binary),str(dest),str(seed),'0.1',backend]
                with (root/f'{name}.stdout').open('w') as out,(root/f'{name}.stderr').open('w') as err:
                    subprocess.run(cmd,env=env,stdout=out,stderr=err,check=True,timeout=660)
                row=json.loads((root/f'{name}.stdout').read_text());row['thread_control']=audit_log((root/f'{name}.stderr').read_text())
                assert row['thread_control']['calls']>=row['finufft_calls']
                if backend=='finufft':assert row['finufft_calls']==row['rho_calls'] and row['thread_control']['calls']>0
                else:assert row['thread_control']['calls']==0
                previous=old/f'five_4t_{seed}_{backend}_{repeat}'
                assert (dest/'initial.tsv').read_bytes()==(previous/'initial.tsv').read_bytes()
                assert (dest/'modes.tsv').read_bytes()==(previous/'modes.tsv').read_bytes()
                row.update(repeat=repeat,initial_sha256=digest(dest/'initial.tsv'),modes_sha256=digest(dest/'modes.tsv'),command=cmd)
                save(root/f'{name}.json',row);rows.append(row)
                print(name,round(row['optimization_seconds'],3),row['direct_reference']['ground_state'],flush=True)
                save(root/'summary.json',dict(status='running',runs=rows,comparisons=pairs))
        for backend in ['direct','finufft']:
            comparison=compare(root/f'{seed}_{backend}_0/final.tsv',root/f'{seed}_{backend}_1/final.tsv',2000)
            pairs.append(dict(seed=seed,comparison=backend+'_repeats',**comparison))
        for repeat in range(2):
            pairs.append(dict(seed=seed,repeat=repeat,comparison='cross_backend',
                **compare(root/f'{seed}_direct_{repeat}/final.tsv',root/f'{seed}_finufft_{repeat}/final.tsv',2000)))
    accepted=all(r['direct_reference']['ground_state'] for r in rows)
    pass_codes=all(0<p['status']<5 for r in rows for p in r['passes'])
    same={b:all(p['bitwise_identical'] for p in pairs if p['comparison']==b+'_repeats') for b in ['direct','finufft']}
    save(root/'summary.json',dict(status='complete',all_ground_states=accepted,all_optimizer_codes_success=pass_codes,
         repeat_endpoints_bitwise_identical=same,runs=rows,comparisons=pairs))
    lines=['# Five-component FINUFFT internal-thread control','',
        'N=2000, five equal components, chi_i=80/798, phi_total=0.2. Fixed-order contact sum. '
        'OMP_NUM_THREADS=OMP_THREAD_LIMIT=4; CPU affinity 0,2,4,6. Only the options passed '
        'to finufft2d1 or finufft_makeplan are changed: opts.nthreads=1. The link wrapper records every override '
        'and status. Original object files, optimizer, tolerance 1e-12, and other options are unchanged.', '',
        'Three seeds, two separate process repeats per backend. Initial files and mode files '
        'match the preceding four-thread FINUFFT experiment byte for byte. Times exclude '
        'independent endpoint validation. Displacements are Cartesian, fixed-label, periodic, '
        'and translation-aligned.', '',
        '| Seed | Vanilla seconds r0/r1 | FINUFFT seconds r0/r1 | Vanilla repeat mean | FINUFFT repeat mean | Cross-backend mean r0/r1 |',
        '|---|---|---|---:|---:|---|']
    for seed in [20260928,20260929,20260930]:
        def times(b):return '/'.join(f"{r['optimization_seconds']:.3f}" for r in sorted(rows,key=lambda r:r['repeat']) if r['seed']==seed and r['backend']==b)
        group=[p for p in pairs if p['seed']==seed]
        means={p['comparison']:p['mean'] for p in group if p['comparison']!='cross_backend'}
        cross='/'.join(f"{p['mean']:.6e}" for p in group if p['comparison']=='cross_backend')
        lines.append(f"| {seed} | {times('direct')} | {times('finufft')} | {means['direct_repeats']:.6e} | {means['finufft_repeats']:.6e} | {cross} |")
    lines += ['',f'All ground-state checks passed: {accepted}. All optimizer codes successful: {pass_codes}.',
              f'Repeat endpoints byte-identical: {same}.', '',
              '| Backend | Median optimization seconds | Minimum final energy | Maximum final energy |',
              '|---|---:|---:|---:|']
    for backend in ['direct','finufft']:
        subset=[r for r in rows if r['backend']==backend];energies=[r['direct_reference']['total_energy'] for r in subset]
        lines.append(f"| {backend} | {statistics.median(r['optimization_seconds'] for r in subset):.6f} | {min(energies):.6e} | {max(energies):.6e} |")
    lines += ['', 'Same-backend reproducibility and agreement between backends are distinct tests. '
              'This finite control does not trace the amplification mechanism or guarantee repeatability '
              'on other machines or configurations. No production thread default was changed; the '
              'broader sweep remains stopped.','']
    (root/'REPORT.md').write_text('\n'.join(lines))
    print('RESULT',json.dumps(dict(all_ground_states=accepted,repeat_identical=same)),flush=True)


if __name__=='__main__':main()
