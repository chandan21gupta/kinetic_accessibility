# Beyond MFE: A kinetics-aware benchmark for RNA inverse design

**Research proposal and reproducible experimental protocol --- v1.0**

## 1. Executive summary

RNA inverse design seeks a sequence that adopts a specified RNA secondary structure. Conventional benchmarks often declare success when a forward folding model predicts a minimum-free-energy (MFE) structure identical to the target. This is a thermodynamic criterion, not a demonstration that the molecule reaches the target during folding. Sequential emergence of RNA can favor alternative structures that persist over experimentally relevant timescales.

We propose to evaluate existing computationally successful designs using cotranscriptional folding simulations. The initial experiment uses synthetic, pseudoknot-free Eterna100-v2 targets and their published solution sequences, ViennaRNA for equilibrium calculations, and DrTransformer for kinetic simulations. The subsequent benchmark compares inverse design algorithms using multiple independent candidates per target. No new inverse design model is required for the pilot.

**Research question:** Among sequences that pass conventional MFE-based inverse design evaluation, how often is the intended target the dominant predicted kinetic structure, and how much population occupies the target over time?

**Hypothesis:** A subset of MFE-successful designs will exhibit low target occupancy at finite times, including cases where the equilibrium target probability is high. Such discrepancies motivate kinetics-aware evaluation, but simulations alone cannot establish experimental folding failure.

## 2. Motivation and literature

RNAInvBench \[1\] evaluates predicted secondary structures against targets using dot-bracket Hamming distance; exact success requires distance zero. Eterna100-v2 \[2\] is a curated set of 100 synthetic design challenges, classified as pseudoknot-free by RNAInvBench (Appendix A.4, Table 4; Appendix A.5). Lengths span 11--399 nt. The revised puzzles were constructed to be solvable under ViennaRNA 2; this does not imply kinetic accessibility. Published reference solutions are reportedly available for 99 targets; verify the actual release and count before analysis.

DrTransformer \[3\] models cotranscriptional folding using a dynamically constructed, coarse-grained secondary-structure energy landscape and population kinetics. Figure 5 of its paper compares the most populated structure with the MFE structure at finite post-transcription times. We adapt that idea by comparing against the original inverse-design target instead of the MFE structure. This is a proposed extension, not a claim that no prior kinetics-aware inverse design studies exist. Conduct a systematic literature review before claiming novelty.

## 3. Definitions and success criteria

For target dot-bracket string T of length L, designed sequence x, and predicted structure S, define HD(S,T) as the number of differing characters. Use canonical pseudoknot-free dot-bracket strings of equal length. HD=0 means exact structural identity. Nonzero HD is not a reliable count of differing base pairs.

-   **Thermodynamic success:** HD(MFE(x),T)=0 under a specified ViennaRNA configuration.
-   **Kinetic dominant-state success at time t:** HD(argmax_S p(S,t),T)=0. Report ties separately; use strict unique dominance as the primary binary rule.
-   **Exact target occupancy:** Pkin(t)=sum of reported occupancies of states whose representative structure equals T.
-   **Off-target representative occupancy:** 1−Pkin(t), provided the distribution is complete and normalized; otherwise report missing mass separately.
-   **Equilibrium target probability:** Peq=exp((Gensemble−Etarget)/(RT)), with energies in kcal/mol, R=0.001987204258 kcal/mol/K and absolute temperature in kelvin.
-   **Equilibrium--kinetic gap:** Δ(t)=Peq−Pkin(t). Positive values indicate a predicted finite-time occupancy deficit relative to equilibrium, not proof of a trap.
-   **Occupancy-weighted normalized HD:** sum_S p(S,t) HD(S,T)/L.
-   **Occupancy-weighted base-pair F1:** sum_S p(S,t) F1(S,T), where F1=2\|B(S)∩B(T)\|/(\|B(S)\|+\|B(T)\|). If both pair sets are empty, define F1=1.

Exact matching is the primary literature-aligned criterion. Supplementary near-target occupancy may use HD/L ≤ 0.05 or F1 ≥ 0.90 and 0.95; these are proposed sensitivity thresholds, not established biological definitions. Distinguish target presence, dominance, majority occupancy and first arrival. A snapshot does not establish first-passage time or persistent folding.

## 4. Phase I: reference-solution pilot

### 4.1 Environment and inputs

Use a Linux environment with Python, ViennaRNA Python bindings, DrTransformer and Git. Example setup (verify package compatibility and command options for the installed release):

``` bash
python -m venv .venv
source .venv/bin/activate
python -m pip install ViennaRNA drtransformer pandas numpy
python -m pip freeze > requirements-lock.txt
DrTransformer --help > drtransformer-help.txt
git clone https://github.com/eternagame/eterna100-benchmarking.git
git -C eterna100-benchmarking rev-parse HEAD > dataset-commit.txt
```

Consult the repository README to locate the v2 puzzle and solution files; do not assume a particular filename. Record the source path, release/commit, file checksum and license. Convert the available reference solutions into `inputs.tsv` with columns `target_id`, `target`, `sequence`, `method`, `candidate_id`. Set method to `reference`. Preserve missing solutions in an exclusion manifest. Uppercase sequences, convert T to U only if documented, validate A/C/G/U alphabet, balanced parentheses, equal lengths, and absence of pseudoknot notation. Do not silently repair invalid structures. Begin with targets shorter than 100 nt for debugging, then run all available valid reference solutions.

### 4.2 Equilibrium baseline

Use the same temperature and nearest-neighbor energy parameter set in both tools where supported (proposed baseline: 37 °C, Turner 2004, standard ViennaRNA defaults otherwise). Record dangling-end settings, lonely-pair rules, salt settings and any differences. Verify compatibility with the published reference benchmark; report rather than hide mismatches.

Illustrative Python calculation using ViennaRNA bindings:

``` python
import RNA, math
md = RNA.md()
md.temperature = 37.0
fc = RNA.fold_compound(sequence, md)
mfe_structure, mfe_energy = fc.mfe()
fc.exp_params_rescale(mfe_energy)
_, ensemble_energy = fc.pf()
target_energy = fc.eval_structure(target)
peq = math.exp((ensemble_energy-target_energy)/(0.001987204258*310.15))
hd = sum(a != b for a,b in zip(mfe_structure,target))
```

Validate equal lengths beforehand. Cross-check Peq against `fc.pr_structure(target)` if available. Investigate invalid energies or probabilities outside \[0,1\] beyond numerical tolerance. Record MFE ties as a limitation of single-structure prediction. Keep all baseline results, but define the primary kinetic cohort as exact MFE matches; report exclusions and analyze nonmatches separately.

### 4.3 Kinetic simulations

Run DrTransformer once per sequence per condition; it estimates populations deterministically rather than requiring 100 stochastic trajectories. Start with the documented default settings, then explicitly record transcription speed, temperature, post-transcription relaxation duration, pruning/coarse-graining parameters and software version. A proposed baseline is 25 nt/s (0.04 s per nucleotide), with observations at completion and 0.04 s, 1 s, 60 s and 3600 s after completion. These are modeling choices, not universal experimental conditions.

Example smoke test (check exact CLI and output filenames against installed help and official documentation):

``` bash
mkdir -p runs/smoke
cd runs/smoke
printf '>test\nGGGGAAAACCCC\n' | DrTransformer > stdout.log 2> stderr.log
```

For production, configure the installed release to retain detailed time-dependent output (.drf or documented equivalent) and simulate at least 3600 s after completion. If the output does not include the requested observation times, use supported time-grid options or separate runs with the required final relaxation durations; do not silently interpolate across transcription events. Store each candidate/condition in a separate directory and save command, raw output, exit status, runtime and settings. Use time relative to full transcript completion, not absolute simulation time. Only compare full-length structures to the full-length target. Prefix states can be inspected separately but are not counted as full-target failures.

### 4.4 Parsing and validation

The documented log contains representative structures, energies, occupancies and IDs; detailed trajectory output provides time-dependent populations. Inspect the actual installed version before implementing a parser. Convert output into a long table with `target_id,candidate_id,condition,time_after_completion_s,transcript_length,structure_id,structure,energy,occupancy`. Verify units and whether occupancy is a fraction or percentage. Check that each structure is valid and has length L at each full-length snapshot. Sum occupancies and report any residual mass; do not silently renormalize incomplete outputs. If a target representative is absent, report zero *reported* target occupancy, not proof of zero physical probability. Distinguish absent target from missing or failed simulations.

For each valid snapshot, compute the metrics in Section 3. Example core logic for a complete normalized distribution:

``` python
assert all(len(s)==len(target) for s,p in states)
assert all(0 <= p <= 1 for s,p in states)
assert abs(sum(p for s,p in states)-1) < 1e-4
hd = lambda s: sum(a != b for a,b in zip(s,target))
pkin = sum(p for s,p in states if hd(s)==0)
maxp = max(p for s,p in states)
dominant = [s for s,p in states if abs(p-maxp)<1e-8]
success = len(dominant)==1 and dominant[0]==target
weighted_hd = sum(p*hd(s)/len(target) for s,p in states)
gap = peq-pkin
```

Aggregate duplicate representative strings before computing dominance. The tolerance above is illustrative; adapt to documented output precision and record it. Validate parser against a manually inspected small run and synthetic test cases: exact target only, target absent, mixed population, tied dominant states, incomplete mass and failed run.

## 5. Controls and interpretation

  ----------------------------------------------------------------------------------------------------------------------------
  Equilibrium target probability   Kinetic target occupancy   Interpretation
  -------------------------------- -------------------------- ----------------------------------------------------------------
  High                             High                       Target thermodynamically favored and accessible in simulation

  Low                              Low                        Thermodynamic competition may explain failure

  High                             Low                        Candidate kinetic limitation; inspect time course and barriers

  Low                              High                       Possible transient kinetic enrichment of target
  ----------------------------------------------------------------------------------------------------------------------------

For descriptive stratification, predefine high Peq ≥ 0.8 and low Pkin ≤ 0.2; report sensitivity to these illustrative thresholds and continuous values. Persistent deficits and long-lived alternative states support a trapping interpretation more strongly than a single snapshot. Inspect structural histories and, where available, transition/barrier information. Compare against full-length folding from an unfolded initial state using a compatible kinetic simulator/configuration if available; do not assume DrTransformer supports this control natively. This distinguishes effects of sequential emergence from general slow folding. Vary transcription rates (e.g. 10, 25, 50 nt/s) and numerical pruning settings on a subset; assess whether conclusions survive model changes. Use matched thermodynamic settings. Include short, simple synthetic hairpins as smoke-test controls, not as guaranteed positive outcomes.

DrTransformer coarse-grains states: occupancy assigned to a representative may include other conformations, and pruning may omit the target. Thus representative occupancy is not exact microstate probability. Report this explicitly and, for selected short sequences, consider comparison with a more exhaustive or stochastic kinetic method. No experimental validation is included in this protocol.

## 6. Phase II: compare inverse design algorithms

After the pilot, choose three or more pseudoknot-free methods (e.g. RNAinverse, IncaRNAtion, Meta-LEARNA, subject to installation and compatibility). Use a fixed panel of 20 targets initially, stratified by length and difficulty, then expand to all 100. Generate 50--100 candidates per method per target with recorded seeds and a predefined computational budget. Record attempted designs, failures, duplicates, timeouts, runtime and model settings. Do not keep generating until a fixed number of successes without reporting the extra cost. Evaluate all valid candidates thermodynamically; run kinetics on the MFE-successful subset and optionally all candidates for unconditional evaluation. Report unique-sequence analysis and attempt-level analysis separately. Record training overlap where known.

For each method and target, report MFE success fraction, kinetic dominant-target fraction among MFE successes, mean target occupancy, mean equilibrium--kinetic gap, and weighted structural similarity. Also report unconditional joint success (MFE match AND kinetic dominant match divided by attempted designs, with missing simulations explicitly accounted for). Targets with zero MFE successes have undefined conditional kinetic metrics, not zero; report their count. Compare methods on common evaluable targets and show coverage. Macro-average target-level values so easy targets with many successes do not dominate. Use paired target-level bootstrap confidence intervals (e.g. 2000 resamples, fixed seed) and report sample sizes; avoid treating correlated candidates or timepoints as independent. Predefine 60 s as the primary post-transcription endpoint, with other times as secondary analyses.

## 7. Outputs and analysis plan

Suggested directory structure:

``` text
project/
  inputs.tsv
  exclusions.tsv
  requirements-lock.txt
  dataset-commit.txt
  config.json
  runs/<candidate_id>/<condition>/
  results/equilibrium.tsv
  results/states.tsv
  results/metrics.tsv
  results/target_summary.tsv
  results/method_summary.tsv
  figures/
```

Produce: (1) Peq vs Pkin scatter with identity line; (2) target occupancy trajectories; (3) distribution of Δ at each time; (4) MFE vs kinetic success rates by method; (5) target-level heatmap and length-stratified summaries; (6) examples of persistent alternative structures; (7) exclusions and simulation failure report. Include uncertainty and denominators. For Phase I, report descriptive target-level results, not algorithm rankings.

Hypothetical example: five MFE-successful sequences have target occupancies 0.95, 0.20, 0.02, 0.85, 0.00 at 60 s. Mean occupancy is 40.4%; if only the first and fourth have target as dominant state, kinetic success is 40% versus 100% MFE success in this selected cohort. These are illustrative, not measured results.

## 8. Expected contribution and limitations

The benchmark would quantify how often equilibrium-based design success translates into finite-time target accessibility under a common kinetic model. A null result (high kinetic accessibility for most designs) is equally informative. Low occupancy alone is not proof of trapping, and in vitro folding depends on conditions, ionic environment, transcription machinery and interactions not fully represented here. Restrict the study to synthetic nonpathogenic benchmark RNAs. Do not claim biological functionality or experimental validation from simulation. A future experimental study could test selected synthetic designs under matched conditions.

**Implementation status:** This is a research protocol, not a completed experiment or tested software pipeline. The CLI flags for detailed output and time scheduling, actual dataset filenames and reference solution count, parameter compatibility and parser schema must be confirmed against the installed releases. No results have been generated. Proposed rates, thresholds, timepoints and sample sizes are design decisions to preregister or revise before production runs.

## 9. References and resources

1.  Cole et al. (2024). RNAInvBench: Benchmark for the RNA Inverse Design Problem. https://openreview.net/pdf?id=zmzBD74UXj ; code: https://github.com/COLA-Laboratory/RNAInvBench
2.  Koodli et al. Redesigning the Eterna100 for the Vienna 2 folding engine. https://doi.org/10.1101/2021.08.26.457839 ; data/code: https://github.com/eternagame/eterna100-benchmarking
3.  DrTransformer: heuristic cotranscriptional RNA folding using the nearest neighbor energy model (2023). https://doi.org/10.1093/bioinformatics/btad034 ; code: https://github.com/ViennaRNA/drtransformer ; package/output documentation: https://pypi.org/project/drtransformer/
4.  ViennaRNA documentation: https://viennarna.readthedocs.io/
5.  Principles for predicting RNA secondary structure design difficulty (2016). https://pmc.ncbi.nlm.nih.gov/articles/PMC4833017/

References are starting points from the preceding research discussion; verify exact release-specific details before running the production benchmark.
