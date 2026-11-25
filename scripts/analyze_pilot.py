"""Validate full-precision representative populations and summarize the pilot."""
import csv,json,math,collections
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results'; FIG=ROOT/'figures'; FIG.mkdir(exist_ok=True)
TIMES=[0,.04,1,60,3600]; RATES=[10,25,50]; TOL=1e-4; TIE=1e-8

def read(p): return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rows):
    if not rows: return
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
def pairs(s):
    stack=[]; b=set()
    for i,c in enumerate(s):
        if c=='(': stack.append(i)
        elif c==')':
            if not stack: raise ValueError('unbalanced structure')
            b.add((stack.pop(),i))
        elif c!='.': raise ValueError('invalid structure character')
    if stack: raise ValueError('unbalanced structure')
    return b
inputs=read(ROOT/'inputs.tsv'); metrics=[]; audit=[]; alternatives=[]
for inp in inputs:
 target=inp['target']; L=len(target); tb=pairs(target); peq=float(inp['peq'])
 assert 0<=peq<=1+1e-8
 for rate in RATES:
    d=ROOT/'runs'/inp['candidate_id']/f'{rate}nt_s'; status=json.loads((d/'status.json').read_text())
    a=dict(target_id=inp['target_id'],candidate_id=inp['candidate_id'],rate_nt_s=rate,length=L,run_status=status['status'],runtime_s=status['runtime_s'],validation='not_evaluated',snapshots=0,error='')
    if status['status']!='complete': audit.append(a); continue
    local=[]; alts=[]
    try:
        groups=collections.defaultdict(lambda:collections.defaultdict(float))
        for r in read(d/'full_length_states.tsv'):
            t=float(r['time_after_completion_s']); p=float(r['occupancy']); s=r['structure']
            if not math.isfinite(t) or not math.isfinite(p) or not 0<=p<=1: raise ValueError('invalid time/population')
            if len(s)!=L: raise ValueError('non-full-length structure')
            pairs(s); groups[t][s]+=p
        for time,states in sorted(groups.items()):
            mass=sum(states.values())
            if abs(mass-1)>TOL: raise ValueError(f'population mass {mass} at {time}')
            pkin=states.get(target,0); maxp=max(states.values()); dom=[s for s,p in states.items() if abs(p-maxp)<=TIE]
            hd=f1=near90=near95=nearhd=0
            for s,p in states.items():
                sb=pairs(s); h=sum(a!=b for a,b in zip(s,target))/L
                f=2*len(sb&tb)/(len(sb)+len(tb)) if sb or tb else 1
                hd+=p*h; f1+=p*f; near90+=p*(f>=.9); near95+=p*(f>=.95); nearhd+=p*(h<=.05)
            match=next((t for t in TIMES if abs(t-time)<1e-8),None)
            local.append(dict(target_id=inp['target_id'],candidate_id=inp['candidate_id'],rate_nt_s=rate,length=L,time_after_completion_s=time,requested_endpoint=match is not None,peq=peq,pkin=pkin,gap=peq-pkin,target_present=target in states,unique_dominant_target=len(dom)==1 and dom[0]==target,dominant_tie=len(dom)>1,target_majority=pkin>.5,weighted_hd=hd,weighted_bp_f1=f1,near_f1_90=near90,near_f1_95=near95,near_hd05=nearhd,reported_mass=mass,residual_mass=1-mass,n_states=len(states)))
            if match in [60,3600]:
                top=sorted(states.items(),key=lambda x:-x[1])[:3]
                for rank,(s,p) in enumerate(top,1): alts.append(dict(target_id=inp['target_id'],rate_nt_s=rate,time_after_completion_s=match,rank=rank,structure=s,occupancy=p,is_target=s==target,target=target))
        found={float(x['time_after_completion_s']) for x in local if x['requested_endpoint']}
        if any(not any(abs(t-u)<1e-8 for u in found) for t in TIMES): raise ValueError('missing requested endpoint')
        a.update(validation='pass',snapshots=len(local)); metrics.extend(local); alternatives.extend(alts)
    except Exception as e: a.update(validation='failed',error=str(e))
    audit.append(a)
write(OUT/'run_audit.tsv',audit); write(OUT/'trajectories.tsv',metrics)
endpoints=[m for m in metrics if m['requested_endpoint']]; write(OUT/'metrics.tsv',endpoints); write(OUT/'dominant_structures.tsv',alternatives)
valid={r:{a['target_id'] for a in audit if a['rate_nt_s']==r and a['validation']=='pass'} for r in RATES}; common=set.intersection(*valid.values())
rng=np.random.default_rng(20261125)
def ci(values):
 v=np.array(values); boot=v[rng.integers(0,len(v),size=(2000,len(v)))].mean(axis=1)
 return tuple(np.quantile(boot,[.025,.975]))
summ=[]
for cohort in ['available','common']:
 for rate in RATES:
  for t in TIMES:
    ms=[m for m in endpoints if m['rate_nt_s']==rate and abs(m['time_after_completion_s']-t)<1e-8 and (cohort=='available' or m['target_id'] in common)]
    if not ms: continue
    lo,hi=ci([m['pkin'] for m in ms]); n=len(ms)
    summ.append(dict(cohort=cohort,rate_nt_s=rate,time_after_completion_s=t,n_evaluable=n,n_selected=50,mean_pkin=np.mean([m['pkin'] for m in ms]),median_pkin=np.median([m['pkin'] for m in ms]),mean_pkin_ci_low=lo,mean_pkin_ci_high=hi,dominant_successes=sum(m['unique_dominant_target'] for m in ms),dominant_fraction=np.mean([m['unique_dominant_target'] for m in ms]),target_absent=sum(not m['target_present'] for m in ms),mean_peq=np.mean([m['peq'] for m in ms]),mean_gap=np.mean([m['gap'] for m in ms]),high_peq_low_pkin=sum(m['peq']>=.8 and m['pkin']<=.2 for m in ms),mean_weighted_hd=np.mean([m['weighted_hd'] for m in ms]),mean_weighted_bp_f1=np.mean([m['weighted_bp_f1'] for m in ms])))
write(OUT/'rate_summary.tsv',summ)
comparisons=[]
for r1,r2 in [(10,25),(10,50),(25,50)]:
 for t in TIMES:
    lookup={(m['target_id'],m['rate_nt_s']):m for m in endpoints if abs(m['time_after_completion_s']-t)<1e-8}
    diff=[lookup[(i,r2)]['pkin']-lookup[(i,r1)]['pkin'] for i in sorted(common,key=int)]
    if diff:
     lo,hi=ci(diff); comparisons.append(dict(rate_from=r1,rate_to=r2,time_after_completion_s=t,n=len(diff),mean_pkin_difference=np.mean(diff),ci_low=lo,ci_high=hi))
write(OUT/'paired_rate_comparisons.tsv',comparisons)
colors=['#0072B2','#D55E00','#009E73']
fig,axs=plt.subplots(1,3,figsize=(13,4),sharex=True,sharey=True)
for ax,rate,c in zip(axs,RATES,colors):
 ms=[m for m in endpoints if m['rate_nt_s']==rate and m['time_after_completion_s']==60]
 ax.scatter([m['peq'] for m in ms],[m['pkin'] for m in ms],s=28,alpha=.7,color=c); ax.plot([0,1],[0,1],'k--',lw=1); ax.set(title=f'{rate} nt/s (n={len(ms)})',xlabel='Equilibrium target probability',xlim=(-.03,1.03),ylim=(-.03,1.03))
axs[0].set_ylabel('Target representative occupancy at 60 s'); fig.tight_layout(); fig.savefig(FIG/'peq_vs_pkin_60s.png',dpi=180); plt.close(fig)
fig,ax=plt.subplots(figsize=(7,4))
for rate,c in zip(RATES,colors):
 ms=[s for s in summ if s['cohort']=='common' and s['rate_nt_s']==rate]
 ax.plot(range(5),[m['mean_pkin'] for m in ms],'o-',label=f'{rate} nt/s',color=c)
 ax.fill_between(range(5),[m['mean_pkin_ci_low'] for m in ms],[m['mean_pkin_ci_high'] for m in ms],alpha=.12,color=c)
ax.set(xticks=range(5),xticklabels=['0','0.04','1','60','3600'],xlabel='Seconds after full-length emergence (categorical spacing)',ylabel='Mean target representative occupancy',ylim=(0,1)); ax.legend(); ax.set_title(f'Common completed targets (n={len(common)}); 95% target bootstrap CI'); fig.tight_layout(); fig.savefig(FIG/'mean_occupancy.png',dpi=180); plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(10,10),sharey=True)
ids=[x['target_id'] for x in inputs]
for ax,rate in zip(axs,RATES):
 lookup={(m['target_id'],m['time_after_completion_s']):m['pkin'] for m in endpoints if m['rate_nt_s']==rate}
 data=np.array([[lookup.get((i,t),np.nan) for t in TIMES] for i in ids]); cmap=plt.colormaps['viridis'].copy(); cmap.set_bad('#cccccc'); im=ax.imshow(data,vmin=0,vmax=1,aspect='auto',cmap=cmap); ax.set(title=f'{rate} nt/s',xticks=range(5),xticklabels=['0','.04','1','60','3600']); ax.tick_params(axis='x',labelrotation=90)
axs[0].set(yticks=range(50),yticklabels=ids,ylabel='Eterna V2 puzzle ID'); fig.colorbar(im,ax=axs.ravel().tolist(),label='Target representative occupancy (gray = unavailable)',fraction=.03); fig.savefig(FIG/'target_heatmap.png',dpi=180,bbox_inches='tight'); plt.close(fig)
print('VALIDATION',collections.Counter(a['validation'] for a in audit),'common targets',len(common))
for s in summ:
 if s['time_after_completion_s']==60: print(s)
print('PAIRED 60s',[x for x in comparisons if x['time_after_completion_s']==60])
print('TIMEOUTS',[(a['target_id'],a['rate_nt_s']) for a in audit if a['run_status']=='timeout'])
