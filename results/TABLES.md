# KRAS G12C QSP — results tables

_Generated 2026-11-04 23:09 UTC._

## Table 1. Calibration: model vs published anchors

| target | arm | metric | observed | predicted | z-residual | source |
| --- | --- | --- | --- | --- | --- | --- |
| soto_ORR | sotorasib | ORR | 37.1 | 42.25 |  1.03 | CodeBreaK100 |
| soto_PFS | sotorasib | PFS |  6.8 |  6.49 | -0.31 | CodeBreaK100 |
| adag_ORR | adagrasib | ORR | 43.0 | 45.00 |  0.33 | KRYSTAL-1 |
| adag_PFS | adagrasib | PFS |  6.5 |  7.20 |  0.58 | KRYSTAL-1 |
| doce_ORR | docetaxel | ORR | 13.2 | 24.00 |  1.80 | CodeBreaK200 |
| doce_PFS | docetaxel | PFS |  4.5 |  3.64 | -0.86 | CodeBreaK200 |
| soto_occ | sotorasib | occupancy |  0.9 |  0.88 | -0.35 | translational |

## Table 2. Fitted parameters

| parameter | estimate | lower | upper |
| --- | --- | --- | --- |
| fb |  2.600 | 0.200 |  4.00 |
| kprog |  0.007 | 0.002 |  0.08 |
| kcov_soto | 14.825 | 2.000 | 40.00 |
| kcov_adag | 22.000 | 2.000 | 60.00 |
| kkill |  0.032 | 0.015 |  0.12 |
| dep |  0.330 | 0.150 |  0.95 |
| bypass |  0.850 | 0.500 |  0.99 |
| fr0 |  0.300 | 0.020 |  0.30 |
| kres |  0.004 | 0.002 |  0.02 |
| kkill_doce |  0.009 | 0.002 |  0.06 |
| kprog_doce |  0.016 | 0.004 |  0.05 |

## Table 3. Virtual-population outcomes (ORR %, median PFS months)

| arm | n | ORR | ORR_lo | ORR_hi | medPFS | mean_occ84 |
| --- | --- | --- | --- | --- | --- | --- |
| sotorasib | 1000 | 41.6 | 38.6 | 44.5 | 6.3 | 0.88 |
| adagrasib | 1000 | 44.4 | 41.3 | 47.5 | 6.9 | 0.91 |
| docetaxel | 1000 | 24.5 | 22.1 | 27.1 | 3.2 | 0.00 |

## Table 4. In silico head-to-head (median [95% CrI])

| comparison | HR | ORR_diff | medPFS_A | medPFS_B |
| --- | --- | --- | --- | --- |
| sotorasib vs docetaxel | 0.69 [0.56, 0.85] | 16.2 [8.7, 23.8] | 6.1 | 3.2 |
| adagrasib vs docetaxel | 0.68 [0.55, 0.83] | 18.3 [11.7, 26.2] | 6.1 | 3.2 |
| sotorasib vs adagrasib | 1.01 [0.83, 1.25] | -2.0 [-10.2, 6.0] | 6.0 | 6.0 |

## Table 5. Practical identifiability (prior->posterior shrinkage)

| parameter | estimate | se | shrinkage | identified |
| --- | --- | --- | --- | --- |
| kcov_adag | 22.000 | 0.1817 | 0.99 | yes |
| kcov_soto | 14.825 | 0.1815 | 0.98 | yes |
| fb |  2.600 | 0.1816 | 0.81 | yes |
| fr0 |  0.300 | 0.0176 | 0.75 | yes |
| dep |  0.330 | 0.1529 | 0.25 | partial |
| bypass |  0.850 | 0.1009 | 0.19 | weak |
| kprog |  0.007 | 0.1556 | 0.00 | no |
| kkill |  0.032 | 0.1808 | 0.00 | no |
| kres |  0.004 | 0.0947 | 0.00 | no |
| kkill_doce |  0.009 | 0.0987 | 0.00 | no |
| kprog_doce |  0.016 | 0.1527 | 0.00 | no |

## Table 6. Global sensitivity (PRCC)

| parameter | soto_ORR | soto_PFS | doce_PFS | signif_band |
| --- | --- | --- | --- | --- |
| kprol | -0.948 | -0.838 | -0.915 | 0.135 |
| kdeath0 |  0.939 |  0.809 |  0.874 | 0.135 |
| dep |  0.522 |  0.253 |  0.056 | 0.135 |
| bypass | -0.410 | -0.175 | -0.071 | 0.135 |
| kkill |  0.411 |  0.141 |  0.024 | 0.135 |
| kfb | -0.067 |  0.110 |  0.059 | 0.135 |
| kcov_soto |  0.156 |  0.108 | -0.018 | 0.135 |
| Kd_sig |  0.122 |  0.061 |  0.066 | 0.135 |
| fb | -0.092 |  0.017 |  0.129 | 0.135 |
| kres | -0.038 | -0.011 | -0.017 | 0.135 |
| fr0 | -0.049 | -0.004 | -0.010 | 0.135 |

## Table 7. Out-of-sample validation

| scheme | held_out | target | observed | predicted | abs_err | pass |
| --- | --- | --- | --- | --- | --- | --- |
| LODO | sotorasib | soto_ORR | 37.1 | 34.57 | 2.53 |  TRUE |
| LODO | sotorasib | soto_PFS |  6.8 |  4.31 | 2.49 | FALSE |
| LODO | sotorasib | soto_occ |  0.9 |  0.87 | 0.03 |  TRUE |
| LODO | adagrasib | adag_ORR | 43.0 | 38.86 | 4.14 |  TRUE |
| LODO | adagrasib | adag_PFS |  6.5 |  4.73 | 1.77 | FALSE |
| LODO | docetaxel | doce_ORR | 13.2 | 22.57 | 9.37 |  TRUE |
| LODO | docetaxel | doce_PFS |  4.5 |  3.04 | 1.46 |  TRUE |
| TEMPORAL(fit ORR->predict PFS) | all PFS | soto_PFS |  6.8 |  4.18 | 2.62 | FALSE |
| TEMPORAL(fit ORR->predict PFS) | all PFS | adag_PFS |  6.5 |  4.48 | 2.02 | FALSE |
| TEMPORAL(fit ORR->predict PFS) | all PFS | doce_PFS |  4.5 |  1.28 | 3.22 | FALSE |

## Table 8. Scenario analysis

| scenario | arm | dose_factor | fb_scale | ORR | ORR_lo | ORR_hi | medPFS | nadir_med |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sotorasib 960 QD (reference) | sotorasib | 1.0 | 1.0 | 42.6 | 39.2 | 46.1 | 7.0 | 85.5 |
| adagrasib 600 BID (reference) | adagrasib | 1.0 | 1.0 | 45.5 | 42.1 | 48.8 | 7.6 | 81.0 |
| sotorasib, +50% exposure | sotorasib | 1.5 | 1.0 | 45.6 | 42.2 | 49.0 | 7.6 | 80.9 |
| sotorasib, -30% exposure | sotorasib | 0.7 | 1.0 | 39.9 | 36.6 | 43.2 | 6.2 | 89.6 |
| sotorasib + feedback blocker (fb x0.3) | sotorasib | 1.0 | 0.3 | 47.1 | 43.5 | 50.5 | 7.9 | 78.1 |
| adagrasib + feedback blocker (fb x0.3) | adagrasib | 1.0 | 0.3 | 48.0 | 44.5 | 51.5 | 8.7 | 73.8 |
| sotorasib +feedback +50% exposure | sotorasib | 1.5 | 0.3 | 48.0 | 44.8 | 51.6 | 8.7 | 73.6 |

