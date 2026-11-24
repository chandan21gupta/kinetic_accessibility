"""Run installed DrTransformer with additional full-length output times; kinetics unchanged."""
import numpy as np
import csv, os
from drtransformer.landscape import TrafoLandscape
from drtransformer.drtransformer import main
original=TrafoLandscape.simulate
def simulate(self,snodes,p0,times,*args,**kwargs):
    # Full transcript is the only step whose length equals the complete sequence.
    if self.transcript_length == len(self.sequence):
        times=np.unique(np.concatenate((times,[0,.04,1,60,3600])))
        kwargs['force']=sorted(set(kwargs.get('force',[]) + [0,.04,1,60,3600]))
        with open(os.environ.get('RNA_FULL_STATES','full_length_states.tsv'),'w') as f:
            w=csv.writer(f,delimiter='\t'); w.writerow(['time_after_completion_s','structure_id','structure','energy','occupancy'])
            for t,pt in original(self,snodes,p0,times,*args,**kwargs):
                for node,p in zip(snodes,pt):
                    w.writerow([format(t,'.17g'),self.nodes[node]['identity'],node,self.nodes[node]['energy']/100,format(float(p),'.17g')])
                f.flush()
                yield t,pt
        return
    yield from original(self,snodes,p0,times,*args,**kwargs)
TrafoLandscape.simulate=simulate
if __name__=='__main__': main()
