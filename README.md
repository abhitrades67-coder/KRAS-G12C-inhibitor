# KRAS G12C inhibitor — Quantitative Systems Pharmacology (QSP) model

A mechanistic QSP model, in R, of **covalent KRAS G12C inhibition** in previously treated
KRAS G12C–mutant non–small-cell lung cancer (NSCLC). It connects drug exposure, covalent
engagement of a cycling target, adaptive MAPK feedback, and tumour kinetics into one framework,
calibrates it to the published endpoints of **sotorasib**, **adagrasib** and a **docetaxel**
comparator, validates it out-of-sample, and uses it for in silico trials and dose/schedule/
combination exploration.

See [`PROTOCOL.md`](PROTOCOL.md) for the full research proposal and pre-specified analysis plan,
and [`verification/SOURCE_VERIFICATION.md`](verification/SOURCE_VERIFICATION.md) for the data sources.

## Scientific scope

- **PK** — one-compartment oral kinetics per drug (sotorasib 960 mg QD, adagrasib 600 mg BID),
  docetaxel as a q3w cytotoxic comparator.
- **Target engagement** — KRAS G12C tracked as GDP-bound (drug-accessible), GTP-bound (active) and
  covalently drug-bound (irreversible) pools; capture competes with nucleotide cycling and protein
  turnover, so newly synthesised (drug-naïve) protein is the engine of escape.
- **Signalling + feedback** — a lumped ppERK output with ERK→upstream negative feedback; inhibiting
  ERK relieves the brake and raises the activation input (adaptive reactivation).
- **Tumour + resistance** — a KRAS-dependent sensitive clone and a drug-tolerant resistant clone,
  with a pre-existing resistant fraction and slow acquired conversion.

## Pipeline

```bash
Rscript run_all.R            # full pipeline (calibration is the long step)
Rscript run_all.R 05 06      # only steps 05* and 06*
```

| Step | Script | Output |
|---|---|---|
| 01 | `scripts/01_verify.R` | `results/verification.csv` (mass balance, monotonicity, unit tests) |
| 02 | `scripts/02_calibrate.R` | `results/fit.rds`, `calibration_*` (DEoptim + Nelder-Mead fit) |
| 03 | `scripts/03_identifiability.R` | `results/identifiability*` (shrinkage, correlations, profiles) |
| 05 | `scripts/05_vpop.R` | `results/vpop_*` (per-patient outcomes, trajectories) |
| 05b | `scripts/05b_ensemble.R` | `results/ensemble_*` (parameter-uncertainty 95% region, propagated) |
| 06 | `scripts/06_insilico_trial.R` | `results/insilico_*` (ORR diff, PFS HR with CrI) |
| 04 | `scripts/04_validate.R` | `results/validation.csv` (leave-one-drug-out, temporal) |
| 07 | `scripts/07_sensitivity.R` | `results/sensitivity_prcc.csv` (LHS + PRCC) |
| 07b | `scripts/07b_structural.R` | `results/structural_comparison.csv` (reduced-model AIC comparison) |
| 55 | `scripts/55_benchmarks.R` | `results/benchmark_*` (MBMA external benchmark) |
| 08 | `scripts/08_scenarios.R` | `results/scenarios.csv` (dose/schedule/combination) |
| 09 | `scripts/09_figures.R` | `figures/*.png` |
| 10 | `scripts/10_tables.R` | `results/TABLES.md`, `results/headline.txt` |

Then knit the report:

```bash
Rscript -e 'rmarkdown::render("REPORT.Rmd")'
```

## Code layout

```
R/            config, params (+ sources), model ODEs, arms, metrics, calibration helpers
scripts/      the numbered pipeline steps above
results/      numbers, tables, fit
figures/      figures
PROTOCOL.md   research proposal + analysis plan
REPORT.Rmd    rendered report
```

## Requirements

R (≥4.3) with `deSolve`, `data.table`, `DEoptim`, `lhs`, `survival`, `numDeriv`, `ggplot2`,
`patchwork`, `scales`, `viridis`, `jsonlite`, `rmarkdown`. On Debian/Ubuntu these are the
`r-cran-*` packages.

## Reproducibility

A single seed (`SEED` in `R/config.R`) drives the virtual population and all resampling; the
calibration objective uses common random numbers so it is smooth and repeatable.

## Disclaimer

Research and educational in silico model. Calibrated to **summary-level published endpoints**, not
individual patient data. **Not** a clinical decision tool and not a treatment recommendation for any
individual.
