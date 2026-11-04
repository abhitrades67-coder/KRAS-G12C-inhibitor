# Model inputs and their sources

Published values used by the model (`R/params.R`), each with its source and its use. Calibration
anchors are **summary-level published endpoints**; no individual patient data are used. Values are
entered as concise point estimates with an approximate uncertainty (`TARGETS` in `R/params.R`);
this file is the provenance record.

## Drugs and regimens

| Input | Value used | Source / basis |
|---|---|---|
| Sotorasib (AMG 510) regimen | 960 mg orally once daily | sotorasib US prescribing information |
| Sotorasib terminal half-life | ~5 h | sotorasib clinical pharmacology |
| Sotorasib apparent volume / protein binding | V/F ~230 L; ~89% bound (fu≈0.11) | sotorasib clinical pharmacology |
| Adagrasib (MRTX849) regimen | 600 mg orally twice daily | adagrasib US prescribing information |
| Adagrasib terminal half-life | ~23 h | adagrasib clinical pharmacology |
| Adagrasib apparent volume / protein binding | V/F ~1100 L; highly bound (fu≈0.02) | adagrasib clinical pharmacology |
| Docetaxel regimen | 75 mg/m² IV every 3 weeks | docetaxel prescribing information; CodeBreaK 200 control arm |

PK parameters set the normalised free exposure that drives covalent capture; because the covalent
potency and the absolute exposure scale are only jointly identifiable, the exposure scale is fixed
from PK and the per-drug covalent rate constant carries the engagement efficiency (see
`results/identifiability_notes.txt`).

## Calibration targets (efficacy and engagement)

| Target (`TARGETS$tkey`) | Value | Source |
|---|---|---|
| `soto_ORR` — sotorasib confirmed ORR, 2L+ NSCLC | 37.1% | CodeBreaK 100 (Skoulidis et al., N Engl J Med 2021) |
| `soto_PFS` — sotorasib median PFS | 6.8 months | CodeBreaK 100; cf. CodeBreaK 200 (de Langen et al., Lancet 2023) ~5.6 months |
| `adag_ORR` — adagrasib confirmed ORR, 2L+ NSCLC | ~43% | KRYSTAL-1 (Jänne et al., N Engl J Med 2022) |
| `adag_PFS` — adagrasib median PFS | 6.5 months | KRYSTAL-1 |
| `doce_ORR` — docetaxel confirmed ORR | ~13% | CodeBreaK 200 control arm |
| `doce_PFS` — docetaxel median PFS | 4.5 months | CodeBreaK 200 control arm |
| `soto_occ` — covalent target occupancy at the clinical dose | ~0.90 | translational / preclinical covalent-occupancy reports |

## Natural history (fixed, not fitted)

| Input | Value used | Basis |
|---|---|---|
| Untreated net tumour growth (SLD) | ~0.007 /day (doubling ~100 days) | aggressive 2L+ NSCLC growth kinetics; `kprol−kdeath0` |
| RECIST assessment interval | 6 weeks | CodeBreaK / KRYSTAL assessment schedules |
| Response definition | confirmed PR/CR: SLD ≤70% of baseline at two consecutive assessments | RECIST v1.1 (confirmation) |
| Progression definition | SLD ≥120% of nadir (+ small absolute increase) | RECIST v1.1 |

## Benchmark data (`scripts/55_benchmarks.R`)

A model-based meta-analysis pools the per-programme ORR point estimates above (random-effects, logit
scale). The response-onset Emax-in-time curve is anchored to the canonical ~6-week median time to
first response in these programmes. The benchmark is an independent, non-mechanistic sanity check.

## Notes on use of published values

Only numerical summary endpoints and standard pharmacology constants are used, each tied to its
public source above. No figures, tables, or text are reproduced from the sources. Where a value is
approximate or pooled across analyses, it is entered with a wider uncertainty in `TARGETS`.
