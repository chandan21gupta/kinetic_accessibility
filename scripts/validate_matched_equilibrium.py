"""Cross-check reconstructed stationary probabilities on small live landscapes."""
import csv,io,json,os,sys,tempfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from plot_matched_parity import stationary
from drtransformer.landscape import TrafoLandscape
from drtransformer.drtransformer import main
original=TrafoLandscape.simulate
checks=[]
def checked(self,snodes,p0,times,*args,**kwargs):
 if self.transcript_length==len(self.sequence):
  pi=stationary([self.nodes[n]['energy']/100 for n in snodes],self.RT)
  actual=self.get_equilibrium_occupancies(snodes)
  error=float(np.max(np.abs(pi-actual))); assert error<1e-10,error
  R=np.zeros((len(snodes),len(snodes)))
  for i,ni in enumerate(snodes):
   for j,nj in enumerate(snodes):
    if self.has_cg_edge(ni,nj): R[j,i]=self.cg_edges[(ni,nj)]['weight']
  np.fill_diagonal(R,-R.sum(axis=0))
  residual=float(np.max(np.abs(R@pi))); relative=residual/max(1,float(np.abs(R).max()))
  assert relative<1e-10,relative
  checks.append(dict(length=len(self.sequence),n_states=len(snodes),max_abs_difference_from_drtransformer=error,relative_stationarity_residual=relative))
 yield from original(self,snodes,p0,times,*args,**kwargs)
TrafoLandscape.simulate=checked
with (ROOT/'inputs.tsv').open() as f: first=next(csv.DictReader(f,delimiter='\t'))['sequence']
with tempfile.TemporaryDirectory(prefix='rna_equilibrium_check_') as d:
 for idx,seq in enumerate(['GGGGAAAACCCC',first]):
  for rate in [10,25,50]:
   sys.stdin=io.StringIO(f'>check\n{seq}\n'); sys.argv=['DrTransformer','--name',f'check{idx}_{rate}','--outdir',d,'--logfile','--t-ext',str(1/rate),'--t-end','3600','-T','37','-d','2']
   main(); checks[-1]['rate_nt_s']=rate
(ROOT/'results/matched_equilibrium_crosscheck.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
