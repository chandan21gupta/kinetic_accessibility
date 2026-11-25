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
fig,axs=plt.subplots(1,3,figsize=(12,3.6),sharey=True)
for ax,rate in zip(axs,['10','25','50']):
 a=[m for m in metrics if m['rate_nt_s']==rate and float(m['time_after_completion_s'])==60]; ax.hist([float(m['gap']) for m in a],bins=np.linspace(-1,1,21)); ax.axvline(0,color='k',ls='--'); ax.set(title=f'{rate} nt/s, n={len(a)}',xlabel='Peq − representative occupancy')
axs[0].set_ylabel('Number of targets'); fig.tight_layout(); fig.savefig(root/'figures/gap_distribution_60s.png',dpi=180); plt.close(fig)
lines=['# Eterna100-V2 kinetic accessibility pilot','', '## Cohort and validation','', 'Fifty reference sequences were sampled with seed 20261125 from 99 valid V2/Vienna2 published solutions. All 99 were exact MFE matches under ViennaRNA 2.8.0 (37 °C, Turner 2004, dangles=2). This is a reference-design pilot, not an algorithm comparison.','', '142/150 simulations completed: 50 at 10 nt/s, 48 at 25 nt/s, and 44 at 50 nt/s. Eight reached the six-hour per-sequence timeout. All completed runs passed structure validity, full-length, finite occupancy, normalization (tolerance 1e-4), and requested-time validation. Duplicate structure strings were aggregated before dominance evaluation (tie tolerance 1e-8). No populations were renormalized. All 50 equilibrium probabilities agree with fc.pr_structure within 1e-5.','', '## Primary endpoint: 60 seconds after full-length emergence','', '| Rate (nt/s) | Evaluable / selected | Dominant target | Mean occupancy | Mean occupancy 95% CI | Low occupancy (≤0.2) |','|---|---|---|---|---|---|']
for s in summary:
 if s['cohort']=='available' and float(s['time_after_completion_s'])==60:
  rate=s['rate_nt_s']; a=[m for m in metrics if m['rate_nt_s']==rate and float(m['time_after_completion_s'])==60]
  lines.append(f"| {rate} | {s['n_evaluable']}/50 | {s['dominant_successes']}/{s['n_evaluable']} ({100*float(s['dominant_fraction']):.1f}%) | {100*float(s['mean_pkin']):.1f}% | {100*float(s['mean_pkin_ci_low']):.1f}–{100*float(s['mean_pkin_ci_high']):.1f}% | {sum(float(m['pkin'])<=.2 for m in a)} |")
lines+=['','The target representative is present in every completed 60-second snapshot, sometimes with very small occupancy. Presence alone does not establish accessibility. There are no cases meeting the protocol’s joint high-equilibrium (Peq ≥0.8), low-kinetic (occupancy ≤0.2) thresholds at 60 seconds.','', '## Paired rate comparison','', 'Forty-four targets completed at all three rates. Their mean 60-second occupancies are 62.2%, 56.0%, and 51.7% at 10, 25, and 50 nt/s. Use this common cohort for rate comparisons; timeouts disproportionately remove challenging targets.','', '| Comparison (faster minus slower) | Mean difference (percentage points) | 95% paired bootstrap CI |','|---|---|---|']
for s in read('results/paired_rate_comparisons.tsv'):
 if float(s['time_after_completion_s'])==60:
  lines.append(f"| {s['rate_to']} − {s['rate_from']} nt/s | {100*float(s['mean_pkin_difference']):.1f} | {100*float(s['ci_low']):.1f} to {100*float(s['ci_high']):.1f} |")
lines+=['','Intervals use 2,000 target-level bootstrap resamples (seed 20261125). They describe variation across these benchmark targets and are not evidence of experimental folding failure. Exploratory comparisons are not adjusted for multiple testing.','', '## Interpretation and limitations','', '- Slower transcription is associated with higher representative target occupancy on the shared cohort. However, cg-auto=10 sets t-fast=t-ext/10: 0.01, 0.004, and 0.002 s at the three rates. Coarse-graining therefore changes with transcription rate. A fixed t-fast sensitivity experiment (e.g. 0.002 s at all three rates) is required before attributing differences to rate alone.','- Average target representative occupancy exceeds the exact equilibrium microstate probability. These quantify different things: DrTransformer merges conformations into representatives and explores a restricted landscape. A negative mean gap is not by itself evidence of transient physical enrichment.','- The final full-length landscape remains fixed during post-transcription relaxation; it does not exhaustively discover new states over the following hour.','- At 3600 seconds the target is dominant in 43/50, 40/48, and 39/44 evaluable simulations. Some alternative representatives persist, but persistent low occupancy alone is insufficient to prove trapping.','- Timeouts: target 6 and 80 at 25 nt/s; targets 6, 73, 74, 80, 86, and 97 at 50 nt/s. These remain missing and must not be coded as zero occupancy.','- ViennaRNA 2.8.0 differs from historical benchmark releases; all selected sequences match their targets in the installed version. No experimental validation is included.','', '## Candidate inspection','', 'At 25 nt/s and 60 s, puzzle 39 has equilibrium target probability 0.776 and target representative occupancy 0.011. It falls just below the preregistered high-equilibrium threshold, so should be inspected as a continuous-value discrepancy rather than labeled a threshold-positive trap. Other low-occupancy candidates include 4, 5, 24, 60, 73, and 74. Inspect their 3600-second states and pruning sensitivity before interpretation.','', '## Outputs','', '- `metrics.tsv`: five requested endpoints for every validated run.','- `trajectories.tsv`: all full-length snapshots.','- `run_audit.tsv`: validation, runtime, and missingness.','- `rate_summary.tsv`: available and common-cohort summaries.','- `paired_rate_comparisons.tsv`: paired occupancy differences.','- `length_summary_60s.tsv`: length-stratified coverage and metrics.','- `dominant_structures.tsv`: top three representatives at 60 and 3600 seconds.','- `../figures/`: equilibrium–kinetic scatter, mean occupancy by endpoint, target heatmap, and gap distribution.','', 'Reproduce with `python scripts/analyze_pilot.py` followed by `python scripts/report_pilot.py`.']
lines.insert(2, '> **Correction:** The original ViennaRNA-versus-DrTransformer scatter and gap compare different state spaces. They are descriptive only, not matched kinetic-deficit measures. Use [the matched-network report](matched_parity_report.md) and `matched_parity.tsv` for the corrected comparison.')
(root/'results/analysis_report.md').write_text('\n'.join(lines)+'\n')
print('Equilibrium crosscheck maximum error',max(errors)); print('Report written')
