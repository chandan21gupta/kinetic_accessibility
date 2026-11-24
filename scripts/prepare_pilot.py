import csv, math, random, hashlib, json
from pathlib import Path
import RNA
root=Path(__file__).resolve().parents[1]
source=root/'data/eterna100-benchmarking/data/eterna100_puzzles.tsv'
RNA.params_load_RNA_Turner2004()
md=RNA.md(); md.temperature=37; md.dangles=2
rows=[]; eligible=[]; exclusions=[]
for r in csv.DictReader(source.open(),delimiter='\t'):
    tid=r['Puzzle #']; target=r['Secondary Structure V2']; seq=r['Sample Solution (V2/Vienna2)'].strip().upper()
    reason=''; depth=0
    if not seq or set(seq)-set('ACGU'): reason='missing_or_invalid_sequence'
    elif len(seq)!=len(target): reason='length_mismatch'
    else:
        for c in target:
            if c not in '.()': reason='invalid_structure'; break
            depth+=(c=='(')-(c==')')
            if depth<0: reason='unbalanced_structure'; break
        if depth: reason='unbalanced_structure'
    if reason:
        exclusions.append(dict(target_id=tid,reason=reason)); continue
    fc=RNA.fold_compound(seq,md); mfe,e=fc.mfe(); fc.exp_params_rescale(e); _,g=fc.pf(); et=fc.eval_structure(target)
    peq=math.exp((g-et)/(.001987204258*310.15))
    hd=sum(a!=b for a,b in zip(mfe,target))
    entry=dict(target_id=tid,name=r['Puzzle Name'],target=target,sequence=seq,length=len(seq),method='reference',candidate_id=f'eterna_v2_{int(tid):03d}',mfe_structure=mfe,mfe_energy=e,target_energy=et,ensemble_energy=g,peq=peq,hd=hd)
    rows.append(entry)
    if hd==0: eligible.append(entry)
    else: exclusions.append(dict(target_id=tid,reason='MFE_nonmatch'))
if len(eligible)<50: raise RuntimeError(f'Only {len(eligible)} eligible')
selected=sorted(random.Random(20261125).sample(eligible,50),key=lambda x:int(x['target_id']))
def write(path,data,fields):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(data)
write(root/'results/equilibrium_all.tsv',rows,list(rows[0]))
write(root/'inputs.tsv',selected,list(selected[0]))
write(root/'exclusions.tsv',exclusions,['target_id','reason'])
config=dict(dataset_commit=(root/'data/eterna100-benchmarking-commit.txt').read_text().strip(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),selection='50 uniformly sampled valid exact MFE matches',seed=20261125,n_valid=len(rows),n_eligible=len(eligible),n_selected=50,viennarna=RNA.__version__,temperature_C=37,dangles=2,parameters='Turner2004',rates_nt_s=[10,25,50],observations_after_full_length_s=[0,.04,1,60,3600],completion_clock='(length-1)*t_ext: first nucleotide emerges at time zero',pruning=.05,cg_auto=10,k0=100000,timeout_per_sequence_s=21600)
(root/'config.json').write_text(json.dumps(config,indent=2)+'\n')
print(json.dumps(config,indent=2)); print('Selected length range',min(x['length'] for x in selected),max(x['length'] for x in selected))
