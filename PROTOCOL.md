# Research proposal and modelling protocol

## A quantitative systems pharmacology (QSP) model of covalent KRAS G12C inhibition in non–small-cell lung cancer: target engagement, adaptive MAPK feedback, resistance, and comparative efficacy

**Version 2.0 — finalised plan (supersedes the v1 outline).**
Modelling is performed entirely in R. This document is the pre-specified analysis plan; every
numbered phase maps to a script under `scripts/` and is orchestrated by `run_all.R`.

---

## 1. Background and rationale

KRAS is the most frequently mutated oncogene in human cancer. The G12C variant — a glycine-to-cysteine
substitution at codon 12 — occurs in roughly 13% of lung adenocarcinomas and a few percent of colorectal
and other tumours. Unlike most KRAS alleles, G12C retains appreciable intrinsic and GAP-stimulated GTP
hydrolysis, so the oncoprotein cycles through a transient GDP-bound (inactive) state. That cycling created
the first therapeutic opening in four decades of "undruggable" KRAS: small molecules that form an
**irreversible covalent bond** with the mutant cysteine-12 in the switch-II pocket, but only while the
protein is GDP-bound, trapping KRAS G12C in its inactive conformation.

Two such inhibitors are approved for previously treated KRAS G12C-mutant NSCLC — **sotorasib** (AMG 510,
960 mg once daily) and **adagrasib** (MRTX849, 600 mg twice daily). Both produce objective responses in a
substantial minority of patients, but responses are frequently **incomplete and short-lived**: median
progression-free survival is on the order of half a year, and the phase III CodeBreaK 200 comparison of
sotorasib against docetaxel showed only a modest progression-free advantage. The central pharmacological
puzzle is **why profferd deep target engagement does not translate into deep, durable tumour control**, and
what model-based changes to dose, schedule, or combination might help.

The mechanistic explanation is an **adaptive feedback** one. ERK, the terminal kinase of the RAS–RAF–MEK–ERK
(MAPK) cascade, normally restrains its own upstream activators (receptor tyrosine kinases, SOS1, and the
phosphatase SHP2). When a G12C inhibitor shuts down ERK, that brake is released within hours to days:
upstream RTK/SOS/SHP2 activity rebounds, newly synthesised (drug-naïve) KRAS G12C is driven into the
active GTP-bound state faster than the covalent drug can capture it, and MAPK output partially
**reactivates** despite the drug still being present. On a longer timescale, **acquired resistance** emerges
from secondary KRAS mutations (e.g. Y96D), bypass signalling, and outgrowth of pre-existing insensitive
clones.

A QSP model is the natural tool to connect these scales — drug exposure, covalent occupancy of a cycling
target, feedback-controlled signalling, and tumour kinetics — into a single quantitative framework that can
be calibrated to clinical endpoints and then interrogated *in silico*.

## 2. Objectives

**Primary.** Build, calibrate, and validate a mechanistic QSP model that reproduces the observed objective
response rate (ORR) and progression-free survival (PFS) of sotorasib and adagrasib in previously treated
KRAS G12C-mutant NSCLC, and of a docetaxel comparator, from first principles of covalent target engagement
and adaptive MAPK feedback.

**Secondary.**
1. Quantify the contribution of adaptive feedback reactivation versus acquired/pre-existing resistance to
   the shape of the tumour-response and PFS curves.
2. Emulate the CodeBreaK 200 head-to-head design *in silico* and recover the sotorasib-versus-docetaxel
   hazard ratio and ORR difference, with credible intervals propagated from a virtual population.
3. Use the model to explore whether alternative schedules (e.g. adagrasib-like continuous twice-daily
   exposure vs once-daily) or a hypothetical feedback-blocking combination (SHP2/SOS1-like) is predicted to
   deepen and prolong response.
4. Identify, by global sensitivity analysis, which biological parameters most strongly govern depth and
   durability of response — i.e. where the model says the therapeutic leverage is.

## 3. Model scope and structure

The model is a system of ordinary differential equations with four coupled modules. State variables and
parameters are defined in `R/params.R`; the right-hand side is `R/model.R`.

### 3.1 Pharmacokinetics (PK)
One-compartment oral PK with first-order absorption for each drug, parameterised from label/literature
clearance, volume, absorption rate, and bioavailability:
- **Sotorasib** 960 mg QD, short half-life (~5 h).
- **Adagrasib** 600 mg BID, long half-life (~23 h), higher accumulation.
- **Docetaxel** 75 mg/m² IV q3w, modelled as a cytotoxic exposure driving a kill term (comparator only).
Free plasma concentration drives target engagement; no tissue compartment is identifiable from the clinical
endpoints used, so the effect-site is taken as free plasma (stated limitation).

### 3.2 Target engagement (covalent, cycling target)
KRAS G12C is tracked in three pools: GDP-bound (inactive, **drug-accessible**), GTP-bound (active,
drug-inaccessible), and **covalently drug-bound** (irreversible). Nucleotide cycling (activation rate set by
upstream RTK/SOS/SHP2 input `U`, deactivation by intrinsic+GAP hydrolysis) competes with covalent capture
(second-order in GDP-pool and free drug, first-order trapping). New KRAS G12C is synthesised drug-naïve and
degraded; this synthesis is the engine of feedback escape. Fractional covalent occupancy is a model output
and a calibration anchor.

### 3.3 MAPK signalling with adaptive feedback
Active KRAS-GTP drives a lumped MAPK output (ppERK, quasi-steady-state Hill response). ERK exerts negative
feedback on the upstream activation input `U` (ERK-dependent suppression of RTK/SOS/SHP2), so inhibition of
ERK relieves the brake and raises `U` — the mechanistic source of reactivation. A dimensionless **pathway
output** `S = ppERK/ppERK_baseline` is the proliferative drive seen by the tumour module.

### 3.4 Tumour growth and resistance
Tumour burden is the sum of a **sensitive** and a **resistant** compartment. Net growth of each = (signal-
driven proliferation − death). The G12C drug reduces `S` in sensitive cells (growth inhibition + an
ERK-suppression-dependent death term); resistant cells carry a reduced drug effect (bypass/secondary
mutation) and a small pre-existing fraction that expands under selection. Observed lesion size (sum of
longest diameters, SLD, as % of baseline) is read out on a RECIST-like schedule to compute ORR and PFS.

## 4. Data and calibration targets

No individual patient data are used. Calibration anchors are **published summary endpoints** from the
pivotal trials, entered as concise cited values in `R/params.R` / `verification/SOURCE_VERIFICATION.md`
(not reproduced here in bulk):

- Sotorasib 960 mg QD, 2L+ NSCLC — ORR and median PFS (CodeBreaK 100; CodeBreaK 200 vs docetaxel).
- Adagrasib 600 mg BID, 2L+ NSCLC — ORR and median PFS (KRYSTAL-1).
- Docetaxel 75 mg/m² — ORR and median PFS (CodeBreaK 200 control arm).
- Target engagement / pharmacodynamics — covalent occupancy and ppERK suppression magnitude/time-course
  from preclinical and translational reports, used as soft priors on the engagement module.

Each target is stored with its point value, an approximate uncertainty, and a source tag.

## 5. Analysis plan (maps to `scripts/` and `run_all.R`)

| Phase | Script | Purpose |
|---|---|---|
| 00 | `00_setup.R` | load packages, source modules, set threads/seed |
| 01 | `01_verify.R` | numerical/code verification: mass balance, QSS checks, dose-event correctness, unit tests |
| 02 | `02_calibrate.R` | fit free parameters to the Section 4 anchors (multi-start + global optimiser); save `fit.rds` |
| 03 | `03_identifiability.R` | prior→posterior shrinkage, profile likelihood, confounded-pair detection |
| 04 | `04_validate.R` | out-of-sample: leave-one-drug-out and temporal (early-visit→late-visit) prediction |
| 05 | `05_vpop.R` | virtual population by LHS, accepted against plausible clinical ranges |
| 05b | `05b_ensemble.R` | parameter-uncertainty ensemble (95% confidence region) propagated into predictions |
| 06 | `06_insilico_trial.R` | emulate CodeBreaK-200-style head-to-head; ORR, PFS, HR with credible intervals |
| 07 | `07_sensitivity.R` | global sensitivity (LHS + PRCC) on depth and durability of response |
| 07b | `07b_structural.R` | structural qualification: refit reduced variants, compare by AIC-like score |
| 08 | `08_scenarios.R` | schedule and combination scenarios (QD vs BID; +feedback blocker) |
| 09 | `09_figures.R` | all figures to `figures/` |
| 10 | `10_tables.R` | all tables + headline numbers to `results/` |

## 6. Validation strategy

The model is considered credible only if it **predicts data it was not fit to**:
1. **Leave-one-drug-out** — refit with one drug's endpoints withheld, predict them from PK + shared mechanism.
2. **Temporal** — fit the early portion of the response/PFS time-course, predict the later portion.
3. **External benchmark** — a model-based meta-analytic (MBMA) time-course of ORR across KRAS G12C programs,
   compared against the mechanistic prediction (`scripts/55_benchmarks.R`).
Pre-specified acceptance: held-out ORR within ±10 percentage points and median PFS within ±1.5 months.

Uncertainty is reported from two sources kept distinct: **patient (trial-to-trial)** variability
(the virtual population and the in-silico trial bootstrap) and **calibration** uncertainty (the
parameter ensemble, `05b`). The per-parameter 95% ensemble range is the marginal
profile-likelihood interval; the full Hessian-based correlation and conditional profiles are in `03`.

## 7. Outputs

`results/` (fitted parameters, calibration table, identifiability, sensitivity indices, in silico trial
summary, validation), `figures/` (PK/PD, target engagement, waterfall, PFS curves, tornado/PRCC,
scenarios), and a rendered `REPORT.html` from `REPORT.Rmd`.

## 8. Assumptions and limitations (pre-stated)

- Effect-site = free plasma (no identifiable tissue compartment from these endpoints).
- MAPK cascade is lumped to a single ppERK output at quasi-steady state; individual kinase layers are not
  resolved.
- Resistance is represented as two compartments (sensitive/resistant) plus a pre-existing fraction, not an
  explicit mutational spectrum.
- Calibration uses summary-level published endpoints, not individual patient data; absolute PFS is matched
  in distribution, not patient-by-patient.
- Docetaxel is a phenomenological cytotoxic comparator, not a mechanistic model of microtubule dynamics.
- This is a research/educational in silico model. It is **not** a clinical decision tool and makes no
  treatment recommendation for any individual.
