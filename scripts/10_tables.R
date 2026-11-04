#!/usr/bin/env Rscript
# 10_tables.R — assemble the headline numbers and a human-readable TABLES.md (PROTOCOL Phase 10).
# Usage: Rscript scripts/10_tables.R
source("R/setup.R"); source("R/calibrate.R")
RES <- DIR_RESULTS
rd <- function(f) if (file.exists(file.path(RES, f))) fread(file.path(RES, f)) else NULL
md <- c("# KRAS G12C QSP — results tables", "",
        paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), " UTC._"), "")
tbl <- function(dt) if (is.null(dt) || !nrow(dt)) "_(not available)_" else {
  h <- paste("|", paste(names(dt), collapse = " | "), "|")
  s <- paste("|", paste(rep("---", ncol(dt)), collapse = " | "), "|")
  rows <- apply(dt, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))
  paste(c(h, s, rows), collapse = "\n")
}

# 1. Calibration
cal <- rd("calibration_targets_vs_pred.csv")
if (!is.null(cal)) cal <- cal[, .(target = tkey, arm, metric, observed = value, predicted = round(pred,2),
                                   `z-residual` = round(resid,2), source)]
md <- c(md, "## Table 1. Calibration: model vs published anchors", "", tbl(cal), "")

# 2. Fitted parameters
fitp <- rd("calibration_fit.csv")
if (!is.null(fitp)) fitp[, estimate := round(estimate, 4)]
md <- c(md, "## Table 2. Fitted parameters", "", tbl(fitp), "")

# 3. Virtual-population outcomes
vs <- rd("vpop_summary.csv")
if (!is.null(vs)) vs <- vs[, .(arm, n, ORR = round(ORR,1), ORR_lo = round(ORR_lo,1), ORR_hi = round(ORR_hi,1),
                               medPFS = round(medPFS,1), mean_occ84 = round(mean_occ84,2))]
md <- c(md, "## Table 3. Virtual-population outcomes (ORR %, median PFS months)", "", tbl(vs), "")

# 4. In silico comparisons
ins <- rd("insilico_summary.csv")
if (!is.null(ins)) ins <- ins[, .(comparison,
   HR = sprintf("%.2f [%.2f, %.2f]", HR, HR_lo, HR_hi),
   ORR_diff = sprintf("%.1f [%.1f, %.1f]", ORRdiff, ORRdiff_lo, ORRdiff_hi),
   medPFS_A = round(medPFS_A,1), medPFS_B = round(medPFS_B,1))]
md <- c(md, "## Table 4. In silico head-to-head (median [95% CrI])", "", tbl(ins), "")

# 5. Identifiability
idt <- rd("identifiability.csv")
if (!is.null(idt)) idt <- idt[, .(parameter, estimate = round(estimate,4), se = round(se,4),
                                  shrinkage = round(shrinkage,2), identified)]
md <- c(md, "## Table 5. Practical identifiability (prior->posterior shrinkage)", "", tbl(idt), "")

# 6. Sensitivity (top drivers)
sp <- rd("sensitivity_prcc.csv")
md <- c(md, "## Table 6. Global sensitivity (PRCC)", "", tbl(sp), "")

# 7. Validation
vl <- rd("validation.csv")
if (!is.null(vl)) vl <- vl[, .(scheme, held_out, key, observed, predicted, abs_err, pass)]
md <- c(md, "## Table 7. Out-of-sample validation", "", tbl(vl), "")

# 8. Scenarios
sc <- rd("scenarios.csv")
md <- c(md, "## Table 8. Scenario analysis", "", tbl(sc), "")

writeLines(md, file.path(RES, "TABLES.md"))

# headline text
hl <- c("KRAS G12C QSP — headline numbers", "================================")
if (!is.null(vs)) for (i in seq_len(nrow(vs)))
  hl <- c(hl, sprintf("%-10s ORR %.1f%% (95%% %.1f-%.1f), median PFS %.1f mo, 12-wk occupancy %.2f",
                      vs$arm[i], vs$ORR[i], vs$ORR_lo[i], vs$ORR_hi[i], vs$medPFS[i], vs$mean_occ84[i]))
if (!is.null(ins)) { hl <- c(hl, "", "In silico comparisons:");
  for (i in seq_len(nrow(ins))) hl <- c(hl, sprintf("  %s: HR %s", ins$comparison[i], ins$HR[i])) }
writeLines(hl, file.path(RES, "headline.txt"))
cat(paste(hl, collapse = "\n"), "\n")
cat("\nSaved results/TABLES.md and results/headline.txt\n")
