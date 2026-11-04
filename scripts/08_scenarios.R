#!/usr/bin/env Rscript
# 08_scenarios.R — model-based exploration of dose/schedule and a feedback-blocking
# combination (PROTOCOL Objective 3). Each scenario is evaluated across the virtual
# population with bootstrap confidence intervals. These are in silico hypotheses, not
# clinical recommendations. Usage: Rscript scripts/08_scenarios.R [N]
source("R/setup.R"); source("R/calibrate.R")
a <- as.numeric(commandArgs(trailingOnly = TRUE)); N <- if (length(a) && !is.na(a[1])) a[1] else 800
par <- load_fit(); vp <- gen_vpop(par, N)

# a scenario modifies the drive and/or the parameters, then we summarise ORR/PFS.
run_scen <- function(label, arm, dose_factor = 1, fb_scale = 1) {
  p <- par; p["fb"] <- par[["fb"]] * fb_scale
  vpl <- lapply(vp, function(v) { v$par["fb"] <- v$par["fb"] * fb_scale; v })
  res <- lapply(vpl, function(v) {
    drv <- arm_drive(arm, v$par, pk_factor = v$pk_factor, dose_factor = dose_factor,
                     kcov_factor = v$kcov_factor, doce_factor = v$doce_factor)
    simulate_patient(v$par, drv) })
  res <- res[!vapply(res, is.null, logical(1))]
  resp <- vapply(res, function(r) best_response(r$SLD), numeric(1))
  pfs  <- vapply(res, function(r) pfs_event(r$time, r$SLD), numeric(1))
  ev   <- as.integer(!is.na(pfs)); tt <- ifelse(is.na(pfs), HORIZON_D, pfs)
  med  <- tryCatch(unname(summary(survfit(Surv(tt,ev)~1))$table["median"])/DAYS_PER_MONTH,
                   error=function(e) median(tt)/DAYS_PER_MONTH); if (is.na(med)) med <- median(tt)/DAYS_PER_MONTH
  orr_ci <- 100*quantile(replicate(600, mean(sample(resp, length(resp), TRUE))), c(.025,.975))
  data.table(scenario = label, arm = arm, dose_factor = dose_factor, fb_scale = fb_scale,
             ORR = round(100*mean(resp),1), ORR_lo = round(orr_ci[1],1), ORR_hi = round(orr_ci[2],1),
             medPFS = round(med,1), nadir_med = round(median(vapply(res,function(r) min(r$SLD),numeric(1))),1))
}

scen <- rbindlist(list(
  run_scen("sotorasib 960 QD (reference)",            "sotorasib"),
  run_scen("adagrasib 600 BID (reference)",           "adagrasib"),
  run_scen("sotorasib, +50% exposure",                "sotorasib", dose_factor = 1.5),
  run_scen("sotorasib, -30% exposure",                "sotorasib", dose_factor = 0.7),
  run_scen("sotorasib + feedback blocker (fb x0.3)",  "sotorasib", fb_scale = 0.3),
  run_scen("adagrasib + feedback blocker (fb x0.3)",  "adagrasib", fb_scale = 0.3),
  run_scen("sotorasib +feedback +50% exposure",       "sotorasib", dose_factor = 1.5, fb_scale = 0.3)
))
fwrite(scen, file.path(DIR_RESULTS, "scenarios.csv"))
cat("== Scenario analysis (ORR %, median PFS months) ==\n"); print(scen)
cat("Saved results/scenarios.csv\n")
