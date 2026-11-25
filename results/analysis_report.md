# Eterna100-V2 kinetic accessibility pilot

The primary evidence is the matched parity plot at 60 seconds. Both axes refer to the same target-labeled representative on the same final retained DrTransformer network, separately for each sequence and transcription condition.

- X: stationary target representative population on that fixed final network.
- Y: kinetic target representative population at 60 seconds after full-length emergence.
- Diagonal: kinetic population equals modeled stationary population.
- Below the diagonal: finite-time target population deficit within this network.
- Above the diagonal: transient enrichment relative to this network's stationary distribution.

The target representative must exactly match the target dot-bracket string. Its population is coarse-grained, not an exact microstate probability.


## Research question and principal observation

Does an exact MFE match guarantee finite-time access to the designed secondary structure? All 50 selected reference sequences have the target as their predicted MFE structure under the installed ViennaRNA configuration. Nevertheless, a subset has high modeled equilibrium target representative population and low target representative population 60 seconds after full-length emergence.

This supports a limitation of an MFE-only evaluation objective: thermodynamic design success does not guarantee rapid target formation within this kinetic model. It does not establish that Eterna100-V2 is generally unreliable, that the sequences fail experimentally, or that the observation is novel relative to prior studies.

## Primary figure

![Matched equilibrium–kinetic populations at 60 seconds](../figures/matched_parity_60s.png)

Each point is one reference sequence. Lower-right points identify candidates whose target representative is favored at modeled equilibrium but weakly populated at the finite observation time. Points near the diagonal show agreement with this network's stationary distribution. High equilibrium preference does not generally imply poor accessibility: many high-population points also lie near the diagonal.

## Reconstruction

The full-precision output saved all representative structures and energies. State sets and energies were checked to remain constant throughout full-length relaxation for all 142 completed runs. DrTransformer 1.0 sets paired transition rates using symmetric saddle energies and representative energies:

`k(i→j) = k0 exp[−(saddle_ij−E_i)/RT]`.

Thus its detailed-balance stationary distribution is `pi_i ∝ exp(−E_i/RT)` over the retained representatives. We used the installed implementation's `RT=0.6163207755 kcal/mol` at 37 °C. No full-ensemble partition function or basin entropy correction is introduced. This is the model's stationary representative distribution, not an exhaustive equilibrium basin probability.

Six small live validation runs (12 and 116 nt, each at three rates) matched DrTransformer's `get_equilibrium_occupancies` to maximum absolute difference 1.7e-16. Stationarity residuals were also checked against the live rate matrices. See `matched_equilibrium_crosscheck.json`.

Original long-run transition matrices and connectivity were not saved, so we cannot retrospectively compute matrix residuals for those exact networks. The reconstruction uses the installed rate formula; completed standard DrTransformer simulations themselves invoke the detailed-balance routine. For future runs, save the final rate matrices and explicit stationary vectors as well.


## Cohort and completed simulations

Fifty sequences were sampled with seed 20261125 from 99 available valid Eterna100-V2/Vienna2 reference solutions, all exact MFE matches under ViennaRNA 2.8.0 (37 °C, Turner 2004, dangles=2). Of 150 attempted rate/sequence simulations, 142 completed and passed full-length structure and population validation, including all five requested times and population sums within 1e-4. No populations were renormalized.

Eight simulations reached the six-hour timeout: puzzles 6 and 80 at 25 nt/s; puzzles 6, 73, 74, 80, 86, and 97 at 50 nt/s. Their populations are missing, not zero. Forty-four targets completed at every rate.

## Results at 60 seconds

| Rate (nt/s) | Completed | Mean stationary target population | Mean kinetic target population | Mean deficit |
|---|---:|---:|---:|---:|
| 10 | 50 | 81.0% | 54.8% | 26.2 percentage points |
| 25 | 48 | 75.4% | 51.4% | 24.0 percentage points |
| 50 | 44 | 71.8% | 51.7% | 20.1 percentage points |

17/50, 17/48, and 12/44 targets have a deficit greater than 0.10 at 60 seconds. This 0.10 cutoff is an illustrative descriptive cutoff added for this plot, not a preregistered trapping definition. At 3600 seconds these counts fall to 12/50, 11/48, and 8/44.

These means are not paired rate comparisons: coverage and retained networks differ across rates. The cg-auto threshold also changes with rate, as discussed in the initial report.

The apparent widespread kinetic enrichment in the original scatter disappears under this matched comparison. Differences between the original and corrected x-values reflect the restricted representative state space as well as coarse-graining; they should not be attributed solely to grouping basin microstates.

## Scope and follow-up

The stationary baseline is conditional on the explored and pruned final network. It can be close to one if few competing representatives survive; in a one-state network it is necessarily one. This does not establish a high equilibrium target probability for the physical molecule. Deficits identify candidates for kinetic investigation within the model, not experimental traps.

Follow up with saved rate matrices, fixed coarse-graining across transcription rates, pruning sensitivity, and selected more exhaustive landscape calculations.

## Files

- `../figures/matched_parity_60s.png` and `.pdf`: primary corrected plot.
- `../figures/matched_parity_0s.*` and `matched_parity_3600s.*`: secondary endpoints.
- `matched_parity.tsv`: matched population metrics at all five requested times.
- `final_network_stationary_states.tsv`: stationary probabilities for each retained representative.
- `matched_parity_audit.tsv`: per-network state-set and energy checks.
- `matched_parity_method.json`: formula, RT, and installed-source checksum.

Reproduce with `python scripts/plot_matched_parity.py` using the local raw outputs. Small-network validation: `python scripts/validate_matched_equilibrium.py`.

## Presentation and reproduction

The main report uses only the matched 60-second parity plot. Original unmatched scatter/gap figures and ancillary descriptive figures are archived under `figures/archive/`; no conclusions about kinetic deficits are drawn from them. Matched 0- and 3600-second figures remain supplementary outputs.

From the project environment, run:

```bash
python scripts/analyze_pilot.py
python scripts/plot_matched_parity.py
python scripts/report_pilot.py
```

The report generator copies this matched-method report into `analysis_report.md` so future regeneration preserves the corrected interpretation.
