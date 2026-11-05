#!/usr/bin/env Rscript
# 05b_ensemble.R — parameter-uncertainty propagation (PROTOCOL Phase 5, uncertainty).
# The point-estimate predictions (05/06) propagate only patient (trial-to-trial) variability.
# Here we build an approximate 95% joint-confidence-region ENSEMBLE of parameter sets
# (Latin-hypercube over the prior box, accepted by the calibration objective via the
# likelihood-ratio threshold obj_min + qchisq(0.95, n_params)), then propagate it into the
# per-arm ORR/PFS and the head-to-head hazard ratio. Identified parameters stay tight; weakly
# identified ones span their prior, so the resulting CrIs honestly reflect calibration
# uncertainty. Usage: Rscript scripts/05b_ensemble.R [M N_acc N_pred cap]
source("R/setup.R"); source("R/calibrate.R")
suppressWarnings(suppressMessages({library(lhs); library(survival)}))
a <- as.numeric(commandArgs(trailingOnly = TRUE))
M     <- if (length(a) >= 1 && !is.na(a[1])) a[1] else 1500   # candidate parameter sets
N_ACC <- if (length(a) >= 2 && !is.na(a[2])) a[2] else 120    # pop size for acceptance objective
N_PRD <- if (length(a) >= 3 && !is.na(a[3])) a[3] else 400    # per-member prediction pool
CAP   <- if (length(a) >= 4 && !is.na(a[4])) a[4] else 150    # members carried into prediction
fit <- readRDS(file.path(DIR_RESULTS, "fit.rds")); best <- fit$best
np  <- length(FIT_PARS)

# objective at the optimum, evaluated at the ensemble's N/seed (consistent threshold)
obj_best <- calib_objective(best, N = N_ACC, seed = SEED, cores = QSP_THREADS)
thresh   <- obj_best + qchisq(0.95, df = np)
cat(sprintf("Ensemble: M=%d  obj_min=%.2f  threshold(95%%, df=%d)=%.2f\n", M, obj_best, np, thresh))

# Latin-hypercube over the prior box; always include the point estimate
set.seed(SEED + 3)
X <- randomLHS(M, np); X <- sweep(sweep(X, 2, (FIT_UB - FIT_LB), `*`), 2, FIT_LB, `+`)
colnames(X) <- FIT_PARS; X <- rbind(best, X)
objs <- unlist(parallel::mclapply(seq_len(nrow(X)), function(i) {
  th <- X[i, ]; names(th) <- FIT_PARS; calib_objective(th, N = N_ACC, seed = SEED, cores = 1) },
  mc.cores = QSP_THREADS))
acc <- which(objs <= thresh)
cat(sprintf("  accepted %d / %d candidates within the 95%% region\n", length(acc), nrow(X)))
if (length(acc) > CAP) acc <- acc[order(objs[acc])][seq_len(CAP)]   # keep the best CAP
ens <- X[acc, , drop = FALSE]; ens_obj <- objs[acc]

# parameter ranges across the ensemble (shows which parameters the data pin vs leave open)
prng <- data.table(parameter = FIT_PARS, estimate = as.numeric(best),
                   ens_lo = apply(ens, 2, quantile, 0.025),
                   ens_med = apply(ens, 2, median),
                   ens_hi = apply(ens, 2, quantile, 0.975))
fwrite(cbind(prng), file.path(DIR_RESULTS, "ensemble_params.csv"))

# propagate each ensemble member: per-arm ORR/PFS and the soto/adag-vs-docetaxel Cox HR
km_med <- function(tt, ev) tryCatch(unname(summary(survfit(Surv(tt, ev) ~ 1))$table["median"])/DAYS_PER_MONTH,
                                    error = function(e) median(tt)/DAYS_PER_MONTH)
member <- function(i) {
  par <- assemble_par(setNames(ens[i, ], FIT_PARS))
  out <- list(); hz <- list()
  arms_d <- list()
  for (arm in ARMS) {
    vp  <- gen_vpop(par, N_PRD, seed = SEED + 100 + i)   # fresh patients per member
    res <- simulate_arm(par, arm, vp, cores = 1)
    resp <- vapply(res, function(r) best_response(r$times, r$sld), numeric(1))
    pfs  <- vapply(res, function(r) pfs_event(r$times, r$sld), numeric(1))
    ev <- as.integer(!is.na(pfs)); tt <- ifelse(is.na(pfs), HORIZON_D, pfs)
    out[[arm]] <- c(ORR = 100*mean(resp), medPFS = km_med(tt, ev))
    arms_d[[arm]] <- data.table(tt = tt, ev = ev)
  }
  hr <- function(A, B) tryCatch(exp(coef(coxph(Surv(tt, ev) ~ grp,
        data = rbind(cbind(arms_d[[A]], grp = 1L), cbind(arms_d[[B]], grp = 0L))))[1]),
        error = function(e) NA_real_)
  data.table(member = i,
    soto_ORR = out$sotorasib["ORR"], soto_PFS = out$sotorasib["medPFS"],
    adag_ORR = out$adagrasib["ORR"], adag_PFS = out$adagrasib["medPFS"],
    doce_ORR = out$docetaxel["ORR"], doce_PFS = out$docetaxel["medPFS"],
    HR_soto_doce = hr("sotorasib","docetaxel"), HR_adag_doce = hr("adagrasib","docetaxel"))
}
ep <- rbindlist(parallel::mclapply(seq_len(nrow(ens)), member, mc.cores = QSP_THREADS))
fwrite(ep, file.path(DIR_RESULTS, "ensemble_predictions.csv"))

qci <- function(x) quantile(x, c(.5, .025, .975), na.rm = TRUE)
summ <- rbindlist(lapply(c("soto_ORR","soto_PFS","adag_ORR","adag_PFS","doce_ORR","doce_PFS",
                           "HR_soto_doce","HR_adag_doce"), function(k) {
  q <- qci(ep[[k]]); data.table(quantity = k, median = q[1], lo95 = q[2], hi95 = q[3]) }))
fwrite(summ, file.path(DIR_RESULTS, "ensemble_summary.csv"))

cat("\n== Parameter ensemble (95% region) ==\n")
print(prng[, .(parameter, estimate = round(estimate,4), ens_lo = round(ens_lo,4), ens_hi = round(ens_hi,4))])
cat("\n== Predictions with PARAMETER uncertainty (median [95% CrI]) ==\n")
print(summ[, .(quantity, median = round(median,2), lo95 = round(lo95,2), hi95 = round(hi95,2))])
cat(sprintf("\nEnsemble size carried into prediction: %d\n", nrow(ens)))
cat("Saved results/ensemble_params.csv, ensemble_predictions.csv, ensemble_summary.csv\n")
