import csv,json,os,subprocess,sys,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
root=Path(__file__).resolve().parents[1]; rate=int(sys.argv[1]); workers=int(sys.argv[2])
rows=list(csv.DictReader((root/'inputs.tsv').open(),delimiter='\t'))
def run(r):
    d=root/'runs'/r['candidate_id']/f'{rate}nt_s'; d.mkdir(parents=True,exist_ok=True)
    fa=d/'input.fa'; fa.write_text('>'+r['candidate_id']+'\n'+r['sequence']+'\n')
    cmd=[sys.executable,str(root/'scripts/drtransformer_observations.py'),'--name',r['candidate_id'],'--outdir',str(d),'--logfile','--t-ext',str(1/rate),'--t-end','3600','-T','37','-d','2','--o-prune','0.05','--cg-auto','10','--k0','100000','--fpwm','4','--mfree','6','--delth','10']
    (d/'command.json').write_text(json.dumps(cmd,indent=2))
    meta=dict(target_id=r['target_id'],rate_nt_s=rate,status='running',start=time.time(),pbs_job_id=os.environ.get('PBS_JOBID'))
    (d/'status.json').write_text(json.dumps(meta,indent=2)); env=os.environ.copy(); env['RNA_FULL_STATES']=str(d/'full_length_states.tsv')
    try:
        with fa.open() as inp,(d/'stdout.log').open('w') as out,(d/'stderr.log').open('w') as err:
            p=subprocess.run(cmd,stdin=inp,stdout=out,stderr=err,env=env,timeout=21600)
        meta.update(status='complete' if p.returncode==0 else 'failed',exit_code=p.returncode)
    except subprocess.TimeoutExpired: meta.update(status='timeout',exit_code=None)
    except Exception as e: meta.update(status='failed',error=str(e))
    meta['runtime_s']=time.time()-meta['start']; (d/'status.json').write_text(json.dumps(meta,indent=2))
    print(r['candidate_id'],rate,meta['status'],flush=True)
    return meta
with ThreadPoolExecutor(max_workers=workers) as pool: statuses=list(pool.map(run,rows))
(root/'results'/f'run_status_{rate}nt_s.json').write_text(json.dumps(statuses,indent=2))
sys.exit(any(s['status']!='complete' for s in statuses))
