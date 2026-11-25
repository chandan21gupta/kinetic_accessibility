# Eterna100-V2 kinetic accessibility pilot

> **Correction:** The original ViennaRNA-versus-DrTransformer scatter and gap compare different state spaces. They are descriptive only, not matched kinetic-deficit measures. Use [the matched-network report](matched_parity_report.md) and `matched_parity.tsv` for the corrected comparison.

## Cohort and validation

Fifty reference sequences were sampled with seed 20261125 from 99 valid V2/Vienna2 published solutions. All 99 were exact MFE matches under ViennaRNA 2.8.0 (37 °C, Turner 2004, dangles=2). This is a reference-design pilot, not an algorithm comparison.

142/150 simulations completed: 50 at 10 nt/s, 48 at 25 nt/s, and 44 at 50 nt/s. Eight reached the six-hour per-sequence timeout. All completed runs passed structure validity, full-length, finite occupancy, normalization (tolerance 1e-4), and requested-time validation. Duplicate structure strings were aggregated before dominance evaluation (tie tolerance 1e-8). No populations were renormalized. All 50 equilibrium probabilities agree with fc.pr_structure within 1e-5.

## Primary endpoint: 60 seconds after full-length emergence

| Rate (nt/s) | Evaluable / selected | Dominant target | Mean occupancy | Mean occupancy 95% CI | Low occupancy (≤0.2) |
|---|---|---|---|---|---|
| 10 | 50/50 | 36/50 (72.0%) | 54.8% | 43.7–65.4% | 13 |
| 25 | 48/50 | 36/48 (75.0%) | 51.4% | 41.3–60.8% | 13 |
| 50 | 44/50 | 36/44 (81.8%) | 51.7% | 41.7–61.1% | 9 |

The target representative is present in every completed 60-second snapshot, sometimes with very small occupancy. Presence alone does not establish accessibility. There are no cases meeting the protocol’s joint high-equilibrium (Peq ≥0.8), low-kinetic (occupancy ≤0.2) thresholds at 60 seconds.

## Paired rate comparison

Forty-four targets completed at all three rates. Their mean 60-second occupancies are 62.2%, 56.0%, and 51.7% at 10, 25, and 50 nt/s. Use this common cohort for rate comparisons; timeouts disproportionately remove challenging targets.

| Comparison (faster minus slower) | Mean difference (percentage points) | 95% paired bootstrap CI |
|---|---|---|
| 25 − 10 nt/s | -6.3 | -10.4 to -2.2 |
| 50 − 10 nt/s | -10.5 | -15.6 to -5.7 |
| 50 − 25 nt/s | -4.2 | -7.9 to -0.9 |

Intervals use 2,000 target-level bootstrap resamples (seed 20261125). They describe variation across these benchmark targets and are not evidence of experimental folding failure. Exploratory comparisons are not adjusted for multiple testing.

## Interpretation and limitations

- Slower transcription is associated with higher representative target occupancy on the shared cohort. However, cg-auto=10 sets t-fast=t-ext/10: 0.01, 0.004, and 0.002 s at the three rates. Coarse-graining therefore changes with transcription rate. A fixed t-fast sensitivity experiment (e.g. 0.002 s at all three rates) is required before attributing differences to rate alone.
- Average target representative occupancy exceeds the exact equilibrium microstate probability. These quantify different things: DrTransformer merges conformations into representatives and explores a restricted landscape. A negative mean gap is not by itself evidence of transient physical enrichment.
- The final full-length landscape remains fixed during post-transcription relaxation; it does not exhaustively discover new states over the following hour.
- At 3600 seconds the target is dominant in 43/50, 40/48, and 39/44 evaluable simulations. Some alternative representatives persist, but persistent low occupancy alone is insufficient to prove trapping.
- Timeouts: target 6 and 80 at 25 nt/s; targets 6, 73, 74, 80, 86, and 97 at 50 nt/s. These remain missing and must not be coded as zero occupancy.
- ViennaRNA 2.8.0 differs from historical benchmark releases; all selected sequences match their targets in the installed version. No experimental validation is included.

## Candidate inspection

At 25 nt/s and 60 s, puzzle 39 has equilibrium target probability 0.776 and target representative occupancy 0.011. It falls just below the preregistered high-equilibrium threshold, so should be inspected as a continuous-value discrepancy rather than labeled a threshold-positive trap. Other low-occupancy candidates include 4, 5, 24, 60, 73, and 74. Inspect their 3600-second states and pruning sensitivity before interpretation.

## Outputs

- `metrics.tsv`: five requested endpoints for every validated run.
- `trajectories.tsv`: all full-length snapshots.
- `run_audit.tsv`: validation, runtime, and missingness.
- `rate_summary.tsv`: available and common-cohort summaries.
- `paired_rate_comparisons.tsv`: paired occupancy differences.
- `length_summary_60s.tsv`: length-stratified coverage and metrics.
- `dominant_structures.tsv`: top three representatives at 60 and 3600 seconds.
- `../figures/`: equilibrium–kinetic scatter, mean occupancy by endpoint, target heatmap, and gap distribution.

Reproduce with `python scripts/analyze_pilot.py` followed by `python scripts/report_pilot.py`.
