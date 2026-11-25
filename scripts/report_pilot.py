"""Create pilot report and supplementary plots from validated metric tables."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import RNA
root=Path(__file__).resolve().parents[1]
def read(p): return list(csv.DictReader((root/p).open(),delimiter='\t'))
inputs=read('inputs.tsv'); metrics=read('results/metrics.tsv'); summary=read('results/rate_summary.tsv'); audit=read('results/run_audit.tsv')
RNA.params_load_RNA_Turner2004(); md=RNA.md(); md.temperature=37; md.dangles=2
errors=[]
for i in inputs:
 fc=RNA.fold_compound(i['sequence'],md); _,e=fc.mfe(); fc.exp_params_rescale(e); fc.pf(); errors.append(abs(fc.pr_structure(i['target'])-float(i['peq'])))
assert max(errors)<1e-5,max(errors)
(root/'results/equilibrium_crosscheck.json').write_text(json.dumps(dict(n=50,max_abs_probability_difference=max(errors),method='fc.pr_structure versus energy-derived probability'),indent=2)+'\n')
rows=[]
for rate in ['10','25','50']:
 for label,low,high in [('<100',0,99),('100–199',100,199),('200–299',200,299),('>=300',300,10000)]:
  a=[m for m in metrics if m['rate_nt_s']==rate and float(m['time_after_completion_s'])==60 and low<=int(m['length'])<=high]
  selected=sum(low<=int(i['length'])<=high for i in inputs)
  rows.append(dict(rate_nt_s=rate,length_bin=label,n_selected=selected,n_evaluable=len(a),mean_pkin=np.mean([float(m['pkin']) for m in a]) if a else '',dominant_successes=sum(m['unique_dominant_target']=='True' for m in a)))
with (root/'results/length_summary_60s.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t'); w.writeheader(); w.writerows(rows)
# Use the matched-network report as the single source for the main analysis.
# This avoids reviving the retired mismatched probability/gap interpretation.
matched = root/'results/matched_parity_report.md'
if not matched.exists():
 raise FileNotFoundError('Generate the matched-network analysis before the main report')
text = matched.read_text().replace('# Matched final-network equilibrium versus kinetic populations', '# Eterna100-V2 kinetic accessibility pilot', 1)
(root/'results/analysis_report.md').write_text(text)
print('Equilibrium crosscheck maximum error',max(errors)); print('Matched-network analysis report written')
