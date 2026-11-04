#!/usr/bin/env Rscript
# 02_calibrate.R — fit the free parameters (FIT_PARS) to the published anchors in TARGETS.
# Robust two-stage search: (1) a large Latin-hypercube scan at small N to locate good basins,
# (2) multi-start Nelder-Mead polish of the best basins at larger N. Common random numbers
# (fixed vpop seed) make the objective deterministic and smooth.
# Usage: Rscript scripts/02_calibrate.R [N M K polish_maxit]
source("R/setup.R"); source("R/calibrate.R")
suppressWarnings(suppressMessages({library(lhs); library(jsonlite)}))

a <- as.numeric(commandArgs(trailingOnly = TRUE))
N      <- if (length(a) >= 1 && !is.na(a[1])) a[1] else 160   # polish/eval population
M      <- if (length(a) >= 2 && !is.na(a[2])) a[2] else 700   # LHS scan points
K      <- if (length(a) >= 3 && !is.na(a[3])) a[3] else 4     # basins to polish
MAXIT  <- if (length(a) >= 4 && !is.na(a[4])) a[4] else 160
NFAST  <- 130
np <- length(FIT_PARS)
cat(sprintf("Calibration: N=%d M=%d K=%d maxit=%d threads=%d (%d params)\n", N, M, K, MAXIT, QSP_THREADS, np))

# ---- stage 1: Latin-hypercube scan (small N, parallel over points) ----
set.seed(SEED)
L <- randomLHS(M, np)
X <- sweep(sweep(L, 2, (FIT_UB - FIT_LB), `*`), 2, FIT_LB, `+`); colnames(X) <- FIT_PARS
X <- rbind(PAR0[FIT_PARS], X)                                 # include the default as a seed
cat("== LHS scan ==\n"); t0 <- Sys.time()
objs <- unlist(parallel::mclapply(seq_len(nrow(X)), function(i) {
  th <- X[i, ]; names(th) <- FIT_PARS; calib_objective(th, N = NFAST, cores = 1) }, mc.cores = QSP_THREADS))
ord <- order(objs)
cat(sprintf("  scanned %d points in %.1f min; best objective (N=%d) = %.2f\n",
            nrow(X), as.numeric(difftime(Sys.time(), t0, units = "mins")), NFAST, objs[ord[1]]))

# ---- stage 2: multi-start Nelder-Mead polish of the top-K basins (larger N) ----
fn <- function(x) { x <- pmin(FIT_UB, pmax(FIT_LB, x)); names(x) <- FIT_PARS; calib_objective(x, N = N, cores = QSP_THREADS) }
cat("== Nelder-Mead polish ==\n")
cands <- list()
for (k in seq_len(K)) {
  start <- X[ord[k], ]; names(start) <- FIT_PARS
  op <- optim(start, fn, method = "Nelder-Mead", control = list(maxit = MAXIT, reltol = 1e-4))
  best <- pmin(FIT_UB, pmax(FIT_LB, op$par)); names(best) <- FIT_PARS
  cands[[k]] <- list(best = best, val = op$value)
  cat(sprintf("  start %d: objective %.3f -> %.3f\n", k, objs[ord[k]], op$value))
}
bi   <- which.min(vapply(cands, function(c) c$val, numeric(1)))
best <- cands[[bi]]$best

# ---- finalise & save (independent evaluation population) ----
par_fit <- assemble_par(best)
pred    <- predict_endpoints(par_fit, N = 400, seed = SEED + 1, cores = QSP_THREADS)
cmp <- copy(TARGETS)[, pred := pred[tkey]][, resid := (pred - value) / sd]
fit_tbl <- data.table(parameter = FIT_PARS, estimate = as.numeric(best),
                      lower = FIT_LB[FIT_PARS], upper = FIT_UB[FIT_PARS])
fwrite(fit_tbl, file.path(DIR_RESULTS, "calibration_fit.csv"))
fwrite(cmp,     file.path(DIR_RESULTS, "calibration_targets_vs_pred.csv"))
saveRDS(list(best = best, par_fit = par_fit, pred = pred, cmp = cmp, obj = cands[[bi]]$val, N = N),
        file.path(DIR_RESULTS, "fit.rds"))
writeLines(toJSON(as.list(round(par_fit, 6)), auto_unbox = TRUE, pretty = TRUE),
           file.path(DIR_CALSTART, "params_fitted.json"))

cat(sprintf("\nFinal objective (eval N=400): %.3f\n", sum(cmp$w * cmp$resid^2)))
print(cmp[, .(tkey, metric, obs = value, pred = round(pred, 2), z = round(resid, 2))])
cat("\nFitted parameters:\n"); print(round(best, 4))
cat("Saved results/calibration_fit.csv, calibration_targets_vs_pred.csv, fit.rds, calibration_start/params_fitted.json\n")
