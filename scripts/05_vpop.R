#!/usr/bin/env Rscript
# 05_vpop.R — virtual population and per-patient outcomes for every arm (PROTOCOL Phase 05).
# Produces the data behind the waterfall, spider and Kaplan-Meier figures and the headline
# ORR/PFS with bootstrap confidence intervals. Usage: Rscript scripts/05_vpop.R [N]
source("R/setup.R"); source("R/calibrate.R")
a <- as.numeric(commandArgs(trailingOnly = TRUE)); N <- if (length(a) && !is.na(a[1])) a[1] else 1000
par <- load_fit()
vp  <- gen_vpop(par, N)
cat(sprintf("Virtual population N=%d; arms: %s\n", N, paste(ARMS, collapse = ", ")))

boot_ci <- function(x, f, B = 1000) {
  v <- replicate(B, f(sample(x, length(x), replace = TRUE))); quantile(v, c(.025, .975), na.rm = TRUE)
}

out <- list(); traj <- list(); summ <- list()
for (arm in ARMS) {
  res <- simulate_arm(par, arm, vp, cores = QSP_THREADS)
  n <- length(res)
  resp  <- vapply(res, function(r) best_response(r$sld), numeric(1))
  nadir <- vapply(res, function(r) min(r$sld), numeric(1))
  pfs   <- vapply(res, function(r) pfs_event(r$times, r$sld), numeric(1))
  ev    <- as.integer(!is.na(pfs)); tt <- ifelse(is.na(pfs), HORIZON_D, pfs)
  best_change <- nadir - 100                                   # % change from baseline (<=−30 is PR)
  med <- tryCatch(unname(summary(survfit(Surv(tt, ev) ~ 1))$table["median"]) / DAYS_PER_MONTH,
                  error = function(e) NA_real_); if (is.na(med)) med <- median(tt) / DAYS_PER_MONTH
  orr <- 100 * mean(resp)
  orr_ci <- 100 * boot_ci(resp, mean)
  out[[arm]] <- data.table(arm = arm, patient = seq_len(n), best_change = best_change,
                           response = resp, pfs_days = tt, pfs_event = ev,
                           occ84 = vapply(res, function(r) r$occ[which.min(abs(r$times - 84))], numeric(1)))
  # keep SLD trajectories for a random subset (spider plot)
  idx <- sort(sample(n, min(60, n)))
  traj[[arm]] <- rbindlist(lapply(idx, function(i)
    data.table(arm = arm, patient = i, time = res[[i]]$times, sld = res[[i]]$sld)))
  summ[[arm]] <- data.table(arm = arm, n = n, ORR = orr, ORR_lo = orr_ci[1], ORR_hi = orr_ci[2],
                            medPFS = med, mean_occ84 = mean(out[[arm]]$occ84))
  cat(sprintf("%-10s ORR=%.1f%% [%.1f,%.1f]  medPFS=%.1f mo  occ=%.2f\n",
              arm, orr, orr_ci[1], orr_ci[2], med, summ[[arm]]$mean_occ84))
}
outcomes <- rbindlist(out); trajectories <- rbindlist(traj); summary_tbl <- rbindlist(summ)

# per-patient parameter table (for sensitivity / phenotype analysis)
ptab <- rbindlist(lapply(seq_along(vp), function(i) {
  pp <- vp[[i]]$par[c("dep","kkill","bypass","fr0","fb","kprol","kres")]
  data.table(patient = i, dep = pp[["dep"]], kkill = pp[["kkill"]], bypass = pp[["bypass"]],
             fr0 = pp[["fr0"]], fb = pp[["fb"]], kprol = pp[["kprol"]], kres = pp[["kres"]],
             pk = unname(vp[[i]]$pk_factor), kcov = unname(vp[[i]]$kcov_factor),
             doce = unname(vp[[i]]$doce_factor)) }))

fwrite(outcomes,     file.path(DIR_RESULTS, "vpop_outcomes.csv"))
fwrite(summary_tbl,  file.path(DIR_RESULTS, "vpop_summary.csv"))
fwrite(ptab,         file.path(DIR_RESULTS, "vpop_patients.csv"))
saveRDS(trajectories, file.path(DIR_RESULTS, "vpop_trajectories.rds"))
cat("Saved results/vpop_outcomes.csv, vpop_summary.csv, vpop_patients.csv, vpop_trajectories.rds\n")
