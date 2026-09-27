#!/usr/bin/env python3
"""Finite five-component queue; prepare validates, run executes the gated study."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("paired", HERE.parent/"finufft_type1/run_comparison.py")
paired=importlib.util.module_from_spec(spec)
spec.loader.exec_module(paired)
REPO=paired.REPO
TARGETS=[0.10,0.20,0.30,0.40,0.45,0.50,0.55]


def stamp():
    return datetime.now(timezone.utc).isoformat()


def write_json(path,value):
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")
    temporary.replace(path)


def environment():
    seen=set();cpus=[]
    for cpu in sorted(os.sched_getaffinity(0)):
        topology=Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        key=((topology/"physical_package_id").read_text(),(topology/"core_id").read_text())
        if key not in seen:
            seen.add(key);cpus.append(cpu)
    if len(cpus)<4:
        raise RuntimeError("four physical cores unavailable")
    os.sched_setaffinity(0,cpus[:4])
    env=dict(os.environ,OMP_NUM_THREADS="4",OMP_THREAD_LIMIT="4",OMP_DYNAMIC="false",
             OMP_WAIT_POLICY="PASSIVE",OMP_PROC_BIND="false",OPENBLAS_NUM_THREADS="4",MKL_NUM_THREADS="4")
    return env,cpus[:4]


def report(root,state):
    rows=[json.loads(p.read_text()) for p in sorted((root/"results").glob("*.json"))] if (root/"results").exists() else []
    lines=["# Five-component FINUFFT study", "", f"Updated: {stamp()}", "", f"Status: **{state['status']}**; phase: {state['phase']}.", "",
           "N=2000 (400/species), 2D, equal radii, number density 1, total packing fraction 0.2. "
           "Four OpenMP threads on distinct physical cores; passive waiting. u0=e0=1. "
           "Complete reciprocal shells; chi is per component, M/(2*399).", "",
           "Optimization uses the existing three-pass L-BFGS routine, at most 2000 evaluations/pass and "
           "600 s/configuration. Reported seconds are actual optimization wall times, not physical simulation time. "
           "Each endpoint is independently recomputed using direct density/force sums and all-pair contact checks.", "",
           "Acceptance: total reference energy<=1e-11, max component S<=1e-10, max force component<=1e-7, "
           "max relative overlap<=1e-6. A passing endpoint and an optimizer success code are distinct recorded facts.", "",
           "## Step 2: paired validation timings", "", "| Seed | Vanilla seconds | FINUFFT seconds | Vanilla accepted | FINUFFT accepted |",
           "|---|---:|---:|---|---|"]
    for seed in [20260928,20260929,20260930]:
        pair={r['backend']:r for r in rows if r['stage']==2 and r['seed']==seed}
        def field(backend,key):
            if backend not in pair:return "pending"
            row=pair[backend]
            return f"{row['optimization_seconds']:.3f}" if key=="time" else str(row['direct_reference']['ground_state'])
        lines.append(f"| {seed} | {field('direct','time')} | {field('finufft','time')} | {field('direct','ok')} | {field('finufft','ok')} |")
    if state.get('step2_positions'):
        lines += ["", "Matched-particle differences after periodic wrapping and global translation:", "",
                  "| Seed | Mean displacement | RMS displacement | Maximum displacement |", "|---|---:|---:|---:|"]
        for seed, values in state['step2_positions'].items():
            lines.append(f"| {seed} | {values['mean']:.4e} | {values['rms']:.4e} | {values['maximum']:.4e} |")
    lines += ["", "## Step 3: FINUFFT-only random-start sweep", "",
              "Vanilla optimization: **not run** in Step 3. Direct summation is used for endpoint validation only.", "",
              "| Target chi_i | Actual chi_i | Completed / 10 | Accepted | FINUFFT median seconds | FINUFFT min–max seconds |",
              "|---|---:|---:|---:|---:|---|"]
    for target in TARGETS:
        group=[r for r in rows if r['stage']==3 and r['target_chi']==target]
        if group:
            times=[r['optimization_seconds'] for r in group]
            lines.append(f"| {target:.2f} | {group[0]['actual_chi']:.6f} | {len(group)} | {sum(r['direct_reference']['ground_state'] for r in group)} | {statistics.median(times):.3f} | {min(times):.3f}–{max(times):.3f} |")
        else:lines.append(f"| {target:.2f} | pending | 0 | — | — | — |")
    lines += ["", "## Individual simulations", "",
              "| Stage | Target chi | Seed | Backend | Seconds | CC energy | Contact energy | max S | max overlap | max force | Accepted | Pass codes |",
              "|---|---:|---|---|---:|---:|---:|---:|---:|---:|---|---|"]
    for row in rows:
        v=row['direct_reference']
        lines.append(f"| {row['stage']} | {row['target_chi']:.2f} | {row['seed']} | {row['backend']} | {row['optimization_seconds']:.3f} | {v['cc_energy']:.3e} | {v['contact_energy']:.3e} | {v['max_S']:.3e} | {v['max_relative_overlap']:.3e} | {v['max_force_component']:.3e} | {v['ground_state']} | {','.join(str(p['status']) for p in row['passes'])} |")
    lines += ["", "## Interpretation and scope", "",
              "Failure to find a candidate within this finite random-start budget does not prove infeasibility. "
              "Validation also tests five translated 20x20 periodic lattices at all seven shell selections; "
              "those constructed candidates are separate from random-start success rates. No chi=0.60, "
              "adaptive refinement, added seeds, or continuation stage is included.", "",
              "Raw receipts, per-run initial/final coordinates and species labels, complete mode lists, "
              "direct-reference structure factors, optimizer codes, and backend-use counters are retained. "
              "state.json records the current job and any blocking condition. The service is a finite Python process, not an agent coordinator.", ""]
    tmp=root/"REPORT.md.tmp";tmp.write_text("\n".join(lines));tmp.replace(root/"REPORT.md")
    write_json(root/"summary.json",{"updated":stamp(),"state":state,"runs":rows})


def update(root,state,status,phase,**extra):
    state.update(status=status,phase=phase,updated=stamp(),**extra)
    write_json(root/"state.json",state);report(root,state)


def positions(root,seed):
    def read(backend):
        file=root/"configurations"/f"step2_chi0.10_seed{seed}_{backend}"/"final.tsv"
        return [tuple(map(float,line.split()[2:4])) for line in file.read_text().splitlines()[1:]]
    a,b=read('direct'),read('finufft');n=len(a);shift=[]
    def wrap(x):return x-math.floor(x+0.5)
    for dim in range(2):
        delta=sorted(wrap(y[dim]-x[dim]) for x,y in zip(a,b));total=math.fsum(delta)
        best=(float('inf'),0)
        for cut in range(n):
            mean=(total+cut)/n
            score=math.fsum(wrap(x-mean)**2 for x in delta)
            if score<best[0]:best=(score,wrap(mean))
        shift.append(best[1])
    distances=[math.sqrt(n)*math.hypot(*(wrap(y[d]-x[d]-shift[d]) for d in range(2))) for x,y in zip(a,b)]
    return dict(mean=statistics.mean(distances),rms=math.sqrt(statistics.mean(x*x for x in distances)),maximum=max(distances),translation_fractional=shift)


def execute(root,state,env,stage,target,seed,backend):
    name=f"step{stage}_chi{target:.2f}_seed{seed}_{backend}"
    receipt=root/"results"/(name+".json")
    if receipt.exists():return json.loads(receipt.read_text())
    directory=root/"configurations"/name
    if directory.exists():raise RuntimeError(f"incomplete prior run {name}; review before resuming")
    update(root,state,"running",f"step{stage}",current_run=name)
    cmd=[state['binary'],str(directory),str(seed),str(target),backend]
    start=stamp()
    with (root/"logs"/(name+".stdout")).open("w") as out,(root/"logs"/(name+".stderr")).open("w") as err:
        proc=subprocess.run(cmd,env=env,cwd=root,stdout=out,stderr=err,timeout=660)
    if proc.returncode:raise RuntimeError(f"{name} failed with exit {proc.returncode}; see logs")
    row=json.loads((root/"logs"/(name+".stdout")).read_text())
    row.update(stage=stage,started=start,finished=stamp(),command=cmd,
               initial_sha256=paired.digest(directory/"initial.tsv"),modes_sha256=paired.digest(directory/"modes.tsv"))
    write_json(receipt,row)
    print(f"{name}: seconds={row['optimization_seconds']:.3f} accepted={row['direct_reference']['ground_state']}",flush=True)
    update(root,state,"running",f"step{stage}",last_completed=name,current_run=None)
    return row


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("action",choices=["prepare","run"]);parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args();root=args.output.resolve();env,cpus=environment()
    if args.action=="prepare":
        root.mkdir(parents=True,exist_ok=False)
        for name in ["results","logs","configurations"]:(root/name).mkdir()
        state={"status":"preparing","phase":"build","created":stamp(),"cpu_affinity":cpus,"N":2000,"packing_fraction":0.2,
               "threads":4,"targets":TARGETS,"seeds_per_target":10,"step2_seeds":[20260928,20260929,20260930],"step3_seeds":list(range(20261001,20261011)),
               "thread_environment":{k:v for k,v in env.items() if k.startswith(("OMP_","OPENBLAS_","MKL_"))},
               "baseline_commit":"7f2fcf535e7213f81c4c4e5709ce93099521340b"}
        update(root,state,"preparing","build")
        app=REPO/"apps/CCOptimization"
        binary=paired.build(root,env,HERE/"experiment.cpp",[app/"MultiStealthyPotential.cpp"],[app])
        state.update(binary=str(binary),binary_sha256=paired.digest(binary))
        update(root,state,"preparing","implementation_validation")
        result=subprocess.run([str(binary),"validate"],env=env,cwd=root,capture_output=True,text=True,timeout=180)
        (root/"logs/validation.stdout").write_text(result.stdout);(root/"logs/validation.stderr").write_text(result.stderr)
        if result.returncode:raise RuntimeError("implementation validation failed; no optimizations authorized")
        validation=json.loads(result.stdout);check=validation['passed']
        if not check:raise RuntimeError("validation did not pass")
        write_json(root/"validation.json",validation)
        update(root,state,"prepared","step2_ready",implementation_validation="passed")
        print(json.dumps(validation,indent=2));return
    state=json.loads((root/"state.json").read_text())
    if state['status'] not in ("prepared","running"):
        raise RuntimeError("study is not prepared/runnable; inspect state")
    if paired.digest(state['binary'])!=state['binary_sha256']:raise RuntimeError("binary changed after validation")
    if cpus!=state['cpu_affinity']:raise RuntimeError("CPU affinity differs from prepared study")
    state['launch_commit']=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip()
    state['runner_sha256']=paired.digest(Path(__file__))
    state['pid']=os.getpid()
    try:
        for seed in state['step2_seeds']:
            pair={}
            for backend in (["direct","finufft"] if seed%2==0 else ["finufft","direct"]):
                pair[backend]=execute(root,state,env,2,0.1,seed,backend)
            if pair['direct']['initial_sha256']!=pair['finufft']['initial_sha256'] or pair['direct']['modes_sha256']!=pair['finufft']['modes_sha256']:
                raise RuntimeError("paired inputs differ")
            state.setdefault('step2_positions',{})[str(seed)]=positions(root,seed)
            if not all(r['direct_reference']['ground_state'] for r in pair.values()):
                update(root,state,"blocked","step2_endpoint_gate",reason="Step 2 did not reach accepted endpoints; sweep not launched")
                return
            if pair['finufft']['optimization_seconds']>=600 or any(p['status']==6 for p in pair['finufft']['passes']):
                update(root,state,"blocked","step2_runtime_gate",reason="FINUFFT exceeded the 10-minute budget; N remains unchanged pending user decision")
                return
        update(root,state,"running","step3",step2="passed",current_run=None)
        for target in TARGETS:
            for seed in state['step3_seeds']:execute(root,state,env,3,target,seed,"finufft")
        update(root,state,"complete","step3_complete",current_run=None)
    except BaseException as error:
        update(root,state,"failed",state['phase'],error=str(error));raise


if __name__=="__main__":main()
