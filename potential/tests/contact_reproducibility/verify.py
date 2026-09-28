"""Verify all paired inputs, endpoint checks, and fixed-sum repeatability receipts."""
import argparse
import json
from pathlib import Path


def main():
    """Check the completed finite matrix and save an explicit pass/fail receipt."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args();root=args.output.resolve()
    checks=[]
    def record(name,value):
        checks.append(dict(check=name,passed=bool(value)))
    count=0
    for variant in ['one_2t','one_4t','five_4t']:
        probe=json.loads((root/'fixed'/f'{variant}_probe.json').read_text())
        direct=next(r for r in probe['backends'] if r['backend']=='direct')
        record(variant+' fixed-point direct values identical',direct['energy_range']==direct['force_range']==0)
        record(variant+' finite difference and reference checks',probe['passed'])
        for seed in [20260928,20260929,20260930]:
            initials=[];modes=[]
            for phase in ['baseline','fixed']:
                rows={}
                for backend in ['direct','finufft']:
                    for repeat in range(2):
                        p=root/phase/f'{variant}_{seed}_{backend}_{repeat}'
                        rows[backend,repeat]=r=json.loads(p.with_suffix('.json').read_text())
                        initials.append((p/'initial.tsv').read_bytes());modes.append((p/'modes.tsv').read_bytes())
                        record(p.name+' '+phase+' independent ground-state check',r['direct_reference']['ground_state'])
                        record(p.name+' '+phase+' optimizer codes',all(0<v['status']<5 for v in r['passes']))
                        record(p.name+' '+phase+' backend used',r['finufft_calls']==(r['rho_calls'] if backend=='finufft' else 0))
                        count+=1
                if phase=='fixed':
                    a=root/phase/f'{variant}_{seed}_direct_0/final.tsv'
                    b=root/phase/f'{variant}_{seed}_direct_1/final.tsv'
                    record(f'{variant} {seed} direct endpoints byte-identical',a.read_bytes()==b.read_bytes())
                    record(f'{variant} {seed} direct pass receipts identical',rows['direct',0]['passes']==rows['direct',1]['passes'])
            record(f'{variant} {seed} all inputs identical across phases/backends/repeats',len(set(initials))==len(set(modes))==1)
    record('72 completed optimizations',count==72)
    receipt=dict(passed=all(x['passed'] for x in checks),optimization_count=count,checks=checks)
    (root/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='checks'}))
    for check in checks:
        if not check['passed']:print('FAILED:',check['check'])
    return 0 if receipt['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
