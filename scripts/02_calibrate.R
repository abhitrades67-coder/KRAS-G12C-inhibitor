#!/usr/bin/env Rscript
# 02_calibrate.R — fit the free parameters (FIT_PARS) to the published anchors in TARGETS.
# Global search (DEoptim) then local polish (Nelder-Mead). Common random numbers (fixed
# vpop seed) make the objective smooth. Usage: Rscript scripts/02_calibrate.R [N NP itermax polish]
source("R/setup.R"); source("R/calibrate.R")
suppressWarnings(suppressMessages({library(DEoptim); library(lhs); library(jsonlite)}))

a <- as.numeric(commandArgs(trailingOnly = TRUE))
N       <- if (length(a) >= 1 && !is.na(a[1])) a[1] else 160
NP      <- if (length(a) >= 2 && !is.na(a[2])) a[2] else 36
ITERMAX <- if (length(a) >= 3 && !is.na(a[3])) a[3] else 22
POLISH  <- if (length(a) >= 4 && !is.na(a[4])) a[4] else 120
cat(sprintf("Calibration: N=%d NP=%d itermax=%d polish=%d threads=%d\n", N, NP, ITERMAX, POLISH, QSP_THREADS))

np <- length(FIT_PARS)
fn <- function(x) { names(x) <- FIT_PARS; calib_objective(x, N = N, cores = QSP_THREADS) }

# seeded initial population: defaults + Latin-hypercube coverage of the box
set.seed(SEED)
L <- randomLHS(NP - 1, np)
ip <- rbind(PAR0[FIT_PARS],
            sweep(sweep(L, 2, (FIT_UB - FIT_LB), `*`), 2, FIT_LB, `+`))
colnames(ip) <- FIT_PARS

cat("== DEoptim global search ==\n"); t0 <- Sys.time()
ctrl <- DEoptim.control(NP = NP, itermax = ITERMAX, strategy = 2, F = 0.7, CR = 0.9,
                        initialpop = ip, trace = 2, parallelType = 0)
de <- DEoptim(fn, lower = FIT_LB, upper = FIT_UB, control = ctrl)
best_de <- de$optim$bestmem; names(best_de) <- FIT_PARS
cat(sprintf("DEoptim best obj=%.3f in %.1f min\n", de$optim$bestval,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cat("== Nelder-Mead polish (N bumped for stability) ==\n")
Npol <- max(N, 220)
fn2 <- function(x) { x <- pmin(FIT_UB, pmax(FIT_LB, x)); names(x) <- FIT_PARS; calib_objective(x, N = Npol, cores = QSP_THREADS) }
op <- optim(best_de, fn2, method = "Nelder-Mead", control = list(maxit = POLISH, reltol = 1e-4))
best <- pmin(FIT_UB, pmax(FIT_LB, op$par)); names(best) <- FIT_PARS

# ---- finalise & save -------------------------------------------------------
par_fit <- assemble_par(best)
pred    <- predict_endpoints(par_fit, N = 400, seed = SEED + 1, cores = QSP_THREADS)  # independent eval pop
cmp <- copy(TARGETS)[, pred := pred[tkey]][, resid := (pred - value) / sd]
fit_tbl <- data.table(parameter = FIT_PARS, estimate = as.numeric(best),
                      lower = FIT_LB[FIT_PARS], upper = FIT_UB[FIT_PARS])

fwrite(fit_tbl, file.path(DIR_RESULTS, "calibration_fit.csv"))
fwrite(cmp,     file.path(DIR_RESULTS, "calibration_targets_vs_pred.csv"))
saveRDS(list(best = best, par_fit = par_fit, pred = pred, cmp = cmp,
             de = de$optim$bestval, obj = op$value, N = N),
        file.path(DIR_RESULTS, "fit.rds"))
writeLines(toJSON(as.list(round(par_fit, 6)), auto_unbox = TRUE, pretty = TRUE),
           file.path(DIR_CALSTART, "params_fitted.json"))

cat(sprintf("\nFinal objective (eval N=400): %.3f\n", sum(cmp$w * cmp$resid^2)))
print(cmp[, .(tkey, metric, obs = value, pred = round(pred, 2), z = round(resid, 2))])
cat("\nFitted parameters:\n"); print(round(best, 4))
cat("Saved results/calibration_fit.csv, calibration_targets_vs_pred.csv, fit.rds, calibration_start/params_fitted.json\n")
