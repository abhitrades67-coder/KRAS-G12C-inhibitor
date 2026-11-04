#!/usr/bin/env Rscript
# 06_insilico_trial.R — emulate a CodeBreaK-200-style randomised head-to-head in silico
# (PROTOCOL Phase 06). Simulates a large pool per arm, then draws B randomised trials to
# propagate trial-to-trial uncertainty into the ORR difference, median PFS and Cox hazard
# ratio (95% credible intervals). Usage: Rscript scripts/06_insilico_trial.R [pool n_per_arm B]
source("R/setup.R"); source("R/calibrate.R")
a <- as.numeric(commandArgs(trailingOnly = TRUE))
POOL <- if (length(a) >= 1 && !is.na(a[1])) a[1] else 4000
NARM <- if (length(a) >= 2 && !is.na(a[2])) a[2] else 300
B    <- if (length(a) >= 3 && !is.na(a[3])) a[3] else 500
par  <- load_fit()
cat(sprintf("In silico trial: pool=%d/arm, n=%d/arm, B=%d trials\n", POOL, NARM, B))

# simulate the outcome pool for each arm
pool <- list()
for (arm in ARMS) {
  vp  <- gen_vpop(par, POOL, seed = SEED + match(arm, ARMS))   # independent pools
  res <- simulate_arm(par, arm, vp, cores = QSP_THREADS)
  pfs <- vapply(res, function(r) pfs_event(r$times, r$sld), numeric(1))
  pool[[arm]] <- data.table(resp = vapply(res, function(r) best_response(r$sld), numeric(1)),
                            pfs_days = ifelse(is.na(pfs), HORIZON_D, pfs),
                            ev = as.integer(!is.na(pfs)))
  cat(sprintf("  pool %-10s ORR=%.1f%% medPFS~%.1f mo\n", arm,
              100*mean(pool[[arm]]$resp), median(pool[[arm]]$pfs_days)/DAYS_PER_MONTH))
}

km_median <- function(d) tryCatch(unname(summary(survfit(Surv(pfs_days, ev) ~ 1, data = d))$table["median"]) / DAYS_PER_MONTH,
                                  error = function(e) median(d$pfs_days)/DAYS_PER_MONTH)

# draw B randomised trials for each pairwise comparison vs docetaxel (+ soto vs adag)
compare <- list(c("sotorasib","docetaxel"), c("adagrasib","docetaxel"), c("sotorasib","adagrasib"))
trials <- list()
set.seed(SEED)
for (cmp in compare) {
  A <- cmp[1]; Bc <- cmp[2]; rows <- vector("list", B)
  for (b in seq_len(B)) {
    da <- pool[[A]][sample(.N, NARM, replace = TRUE)][, grp := 1L]
    db <- pool[[Bc]][sample(.N, NARM, replace = TRUE)][, grp := 0L]
    d  <- rbind(da, db)
    hr <- tryCatch(exp(coef(coxph(Surv(pfs_days, ev) ~ grp, data = d))[1]), error = function(e) NA_real_)
    rows[[b]] <- data.table(comparison = paste(A, "vs", Bc), b = b,
                            ORR_A = 100*mean(da$resp), ORR_B = 100*mean(db$resp),
                            ORR_diff = 100*(mean(da$resp) - mean(db$resp)),
                            medPFS_A = km_median(da), medPFS_B = km_median(db), HR = hr)
  }
  trials[[paste(A, Bc)]] <- rbindlist(rows)
}
trial_dt <- rbindlist(trials)

qci <- function(x) quantile(x, c(.5, .025, .975), na.rm = TRUE)
summ <- trial_dt[, {
  h <- qci(HR); o <- qci(ORR_diff); pa <- qci(medPFS_A); pb <- qci(medPFS_B)
  .(HR = h[1], HR_lo = h[2], HR_hi = h[3],
    ORRdiff = o[1], ORRdiff_lo = o[2], ORRdiff_hi = o[3],
    medPFS_A = pa[1], medPFS_B = pb[1])
}, by = comparison]

fwrite(trial_dt, file.path(DIR_RESULTS, "insilico_trials.csv"))
fwrite(summ,     file.path(DIR_RESULTS, "insilico_summary.csv"))
cat("\n== In silico comparisons (median [95% CrI]) ==\n")
print(summ[, .(comparison, HR = sprintf("%.2f [%.2f,%.2f]", HR, HR_lo, HR_hi),
               ORRdiff = sprintf("%.1f [%.1f,%.1f]", ORRdiff, ORRdiff_lo, ORRdiff_hi),
               medPFS_A = round(medPFS_A,1), medPFS_B = round(medPFS_B,1))])
cat("Saved results/insilico_trials.csv, insilico_summary.csv\n")
