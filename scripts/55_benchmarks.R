#!/usr/bin/env Rscript
# 55_benchmarks.R — external, non-mechanistic benchmark (PROTOCOL Phase 6.3).
# A model-based meta-analysis (MBMA) pooling published objective-response rates across the
# KRAS G12C programmes (random-effects, logit scale), and an Emax-in-time response-onset
# curve. The mechanistic QSP prediction is overlaid and scored (RMSE) against both, as an
# independent sanity check that does not share the QSP's structure.
# Usage: Rscript scripts/55_benchmarks.R
source("R/setup.R"); source("R/calibrate.R")
par <- load_fit()

# --- B1. ORR across programmes (concise, cited point estimates; see SOURCE_VERIFICATION.md) ---
orr_data <- data.table(
  drug  = c("sotorasib","sotorasib","adagrasib"),
  arm   = c("CodeBreaK100","CodeBreaK200","KRYSTAL-1"),
  dose  = c("960 QD","960 QD","600 BID"),
  resp  = c(46, 48, 48),   # approx responders
  n     = c(124, 169, 112),# approx analysed
  orr   = c(37.1, 28.1, 42.9))
orr_data[, `:=`(p = resp/n, logit = qlogis(resp/n), v = 1/(resp) + 1/(n-resp))]
# random-effects pooled ORR (DerSimonian-Laird on logit scale)
w  <- 1/orr_data$v; mu <- sum(w*orr_data$logit)/sum(w)
Q  <- sum(w*(orr_data$logit-mu)^2); tau2 <- max(0,(Q-(nrow(orr_data)-1))/(sum(w)-sum(w^2)/sum(w)))
wr <- 1/(orr_data$v+tau2); mu_r <- sum(wr*orr_data$logit)/sum(wr); se_r <- sqrt(1/sum(wr))
pooled <- plogis(mu_r); pooled_ci <- plogis(mu_r + c(-1,1)*1.96*se_r)
cat(sprintf("MBMA pooled ORR (random-effects) = %.1f%% [%.1f, %.1f]; tau^2(logit)=%.3f\n",
            100*pooled, 100*pooled_ci[1], 100*pooled_ci[2], tau2))

# mechanistic predicted ORR per programme (approved doses)
pred <- predict_endpoints(par, N = 500, seed = SEED + 11, cores = QSP_THREADS)
mech <- data.table(drug = c("sotorasib","adagrasib"),
                   mech_ORR = c(pred[["soto_ORR"]], pred[["adag_ORR"]]))
orr_cmp <- merge(orr_data[, .(drug, arm, obs_ORR = orr)], mech, by = "drug", all.x = TRUE)
orr_cmp[, mbma_pooled := 100*pooled]
rmse_orr <- sqrt(mean((orr_cmp$obs_ORR - orr_cmp$mech_ORR)^2, na.rm = TRUE))
cat(sprintf("Mechanistic vs observed ORR RMSE = %.1f pp\n", rmse_orr))

# --- B2. Response-onset shape: Emax-in-time vs the mechanistic cumulative ORR ---
# canonical onset anchors (cumulative % of eventual responders achieving response by week t)
onset <- data.table(week = c(6, 12, 18, 24), cum_frac = c(0.55, 0.85, 0.95, 1.0))
em <- nls(cum_frac ~ wk^g/(T50^g + wk^g), data = data.table(wk = onset$week, cum_frac = onset$cum_frac),
          start = list(T50 = 7, g = 2), control = nls.control(warnOnly = TRUE))
# mechanistic cumulative ORR over time for sotorasib
vp  <- gen_vpop(par, 400, seed = SEED + 12)
res <- simulate_arm(par, "sotorasib", vp, cores = QSP_THREADS)
wks <- ASSESS_GRID / 7
cumresp <- sapply(seq_along(ASSESS_GRID), function(j) {
  mean(vapply(res, function(r) as.integer(any(r$sld[1:j] <= RECIST_PR)), numeric(1))) })
final <- tail(cumresp, 1); mech_frac <- if (final > 0) cumresp/final else cumresp
onset_cmp <- data.table(week = wks, mbma = predict(em, newdata = data.table(wk = wks)),
                        mech = mech_frac)[week <= 30]

fwrite(orr_cmp,  file.path(DIR_RESULTS, "benchmark_orr.csv"))
fwrite(onset_cmp,file.path(DIR_RESULTS, "benchmark_onset.csv"))
saveRDS(list(pooled = pooled, pooled_ci = pooled_ci, tau2 = tau2, rmse_orr = rmse_orr),
        file.path(DIR_RESULTS, "benchmark_mbma.rds"))
cat("\n== ORR benchmark ==\n"); print(orr_cmp)
cat("Saved results/benchmark_orr.csv, benchmark_onset.csv, benchmark_mbma.rds\n")
