"""Compare target populations on the same final DrTransformer representative network.

DrTransformer 1.0 sets k(i->j)=k0*exp(-(saddle_ij-E_i)/RT), with
symmetric saddle energies. Detailed balance therefore gives stationary weights
pi_i proportional to exp(-E_i/RT) over the retained final representatives.
These are NOT full-ensemble probabilities or sums over basin microstates.
Final energies and all representatives were saved in full_length_states.tsv.
Original transition matrices were not saved; graph residuals cannot be checked
retrospectively. Completed standard simulations use get_p8_detbal internally.
"""
import csv,json,hashlib,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
RT=.61632077549999997 # TrafoLandscape.RT at 37 C, exactly as installed

def read(path):
 with path.open() as f: return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
def stationary(energies,rt=RT):
 e=np.asarray(energies,dtype=float); w=np.exp(-(e-e.min())/rt)
 return w/w.sum()
def main():
 inputs={r['candidate_id']:r for r in read(ROOT/'inputs.tsv')}
 metrics=read(ROOT/'results/metrics.tsv'); output=[]; states_out=[]; audit=[]
 for candidate,inp in inputs.items():
  for rate in (10,25,50):
   d=ROOT/'runs'/candidate/f'{rate}nt_s'; status=json.loads((d/'status.json').read_text())
   if status['status']!='complete': continue
   # At full length, this state set and rate matrix remain fixed over relaxation.
   by_time={}
   with (d/'full_length_states.tsv').open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
     t=float(r['time_after_completion_s']); s=r['structure']; e=float(r['energy']); p=float(r['occupancy'])
     if not all(math.isfinite(x) for x in (t,e,p)): raise ValueError('nonfinite output')
     frame=by_time.setdefault(t,{})
     if s in frame and frame[s][0]!=e: raise ValueError('duplicate structure with conflicting energies')
     old=frame.get(s,(e,0)); frame[s]=(e,old[1]+p)
   initial=by_time[0.0]; names=sorted(initial); energies=[initial[s][0] for s in names]; pi=stationary(energies); eq=dict(zip(names,pi))
   for frame in by_time.values():
    assert set(frame)==set(initial), 'state set changed after full-length emergence'
    assert all(frame[s][0]==initial[s][0] for s in names), 'state energy changed'
    assert abs(sum(p for e,p in frame.values())-1)<1e-4
   target=inp['target']; assert target in eq, 'target not represented'
   for s in names: states_out.append(dict(candidate_id=candidate,target_id=inp['target_id'],rate_nt_s=rate,structure=s,energy_kcal_mol=initial[s][0],stationary_probability=eq[s],is_target=s==target))
   for m in metrics:
    if m['candidate_id']!=candidate or int(m['rate_nt_s'])!=rate: continue
    t=float(m['time_after_completion_s']); kinetic=by_time[t][target][1]
    assert abs(kinetic-float(m['pkin']))<1e-12
    output.append(dict(target_id=inp['target_id'],candidate_id=candidate,rate_nt_s=rate,time_after_completion_s=t,n_final_representatives=len(names),equilibrium_target_representative=eq[target],kinetic_target_representative=kinetic,matched_gap=eq[target]-kinetic,total_variation_from_stationary=.5*sum(abs(by_time[t][s][1]-eq[s]) for s in names),viennarna_exact_target_probability=float(inp['peq'])))
   audit.append(dict(candidate_id=candidate,rate_nt_s=rate,n_states=len(names),state_set_constant=True,energy_constant=True,stationary_mass=float(sum(pi))))
 write(ROOT/'results/matched_parity.tsv',output); write(ROOT/'results/final_network_stationary_states.tsv',states_out); write(ROOT/'results/matched_parity_audit.tsv',audit)
 summary=[]
 for time in (0,60,3600):
  fig,axs=plt.subplots(1,3,figsize=(12,4),sharex=True,sharey=True)
  for ax,rate,color in zip(axs,(10,25,50),('#0072B2','#D55E00','#009E73')):
   rows=[r for r in output if r['rate_nt_s']==rate and r['time_after_completion_s']==time]
   x=np.array([r['equilibrium_target_representative'] for r in rows]); y=np.array([r['kinetic_target_representative'] for r in rows])
   ax.scatter(x,y,s=30,color=color,alpha=.75,edgecolors='white',linewidths=.4)
   ax.plot([0,1],[0,1],'k--',lw=1)
   ax.set(title=f'{rate} nt/s (n={len(rows)})',xlim=(-.03,1.03),ylim=(-.03,1.03),aspect='equal',xlabel='Equilibrium target representative population\n(final retained network)')
   summary.append(dict(rate_nt_s=rate,time_after_completion_s=time,n=len(rows),mean_stationary_target=float(x.mean()),mean_kinetic_target=float(y.mean()),mean_matched_gap=float((x-y).mean()),deficit_over_0p1=int(sum(x-y>.1)),enrichment_over_0p1=int(sum(y-x>.1)),within_0p01=int(sum(abs(x-y)<=.01))))
  axs[0].set_ylabel(f'Kinetic target representative population\nat {time:g} s after full-length emergence')
  fig.suptitle('Same final states and energy model on both axes',fontsize=12)
  fig.tight_layout(); fig.savefig(ROOT/f'figures/matched_parity_{time}s.png',dpi=200); fig.savefig(ROOT/f'figures/matched_parity_{time}s.pdf'); plt.close(fig)
 write(ROOT/'results/matched_parity_summary.tsv',summary)
 source=ROOT/'.conda-env/lib/python3.11/site-packages/drtransformer/landscape.py'
 (ROOT/'results/matched_parity_method.json').write_text(json.dumps(dict(drtransformer='1.0',temperature_C=37,RT_kcal_mol=RT,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),formula='pi_i=exp(-(E_i-min(E))/RT)/sum_j exp(-(E_j-min(E))/RT)',scope='Stationary populations on each run\'s final retained representative network, not full thermodynamic ensemble or explicit basin partition sums',matrix_saved=False,connectivity_verified_retrospectively=False),indent=2)+'\n')
 print('validated networks',len(audit)); print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
