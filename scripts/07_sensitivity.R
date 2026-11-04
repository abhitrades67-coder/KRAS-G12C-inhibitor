#!/usr/bin/env Rscript
# 07_sensitivity.R — global sensitivity of the key clinical outputs to the model parameters
# (PROTOCOL Phase 07). Latin-hypercube sample over plausible parameter ranges, then partial
# rank correlation coefficients (PRCC) of each input against sotorasib ORR, sotorasib median
# PFS and docetaxel median PFS. Usage: Rscript scripts/07_sensitivity.R [n_samples]
source("R/setup.R"); source("R/calibrate.R")
suppressWarnings(suppressMessages(library(lhs)))
a <- as.numeric(commandArgs(trailingOnly = TRUE)); M <- if (length(a) && !is.na(a[1])) a[1] else 220
par0 <- load_fit()

# inputs varied for sensitivity: the fitted mechanism parameters, +/-50% around the fit (within bounds)
sens_pars <- c("fb","kcov_soto","kkill","dep","bypass","fr0","kres","kprol","kdeath0","kfb","Kd_sig")
lo <- sapply(sens_pars, function(p) 0.5 * par0[[p]])
hi <- sapply(sens_pars, function(p) 1.5 * par0[[p]])
# keep the fitted ones inside their calibration bounds
for (p in intersect(sens_pars, FIT_PARS)) { lo[p] <- max(lo[p], FIT_LB[p]); hi[p] <- min(hi[p], FIT_UB[p]) }

set.seed(SEED)
X <- randomLHS(M, length(sens_pars)); colnames(X) <- sens_pars
X <- sweep(sweep(X, 2, (hi - lo), `*`), 2, lo, `+`)

cat(sprintf("Global sensitivity: %d LHS samples x %d inputs\n", M, length(sens_pars)))
eval_one <- function(i) {
  p <- par0; p[sens_pars] <- X[i, ]
  vp <- gen_vpop(p, 90, seed = SEED)                 # small pop; CRN keeps it comparable
  s_soto <- summarise_arm(p, "sotorasib", vp, cores = 1)
  s_doce <- summarise_arm(p, "docetaxel", vp, cores = 1)
  c(soto_ORR = s_soto$ORR, soto_PFS = s_soto$medPFS, doce_PFS = s_doce$medPFS)
}
Y <- do.call(rbind, parallel::mclapply(seq_len(M), eval_one, mc.cores = QSP_THREADS))
ok <- apply(Y, 1, function(r) all(is.finite(r)))
X <- X[ok, ]; Y <- Y[ok, ]
cat(sprintf("  %d/%d samples usable\n", sum(ok), M))

# PRCC: partial correlation of rank(X_j), rank(y) adjusting for the other ranked inputs
prcc <- function(X, y) {
  Rx <- apply(X, 2, rank); ry <- rank(y)
  sapply(seq_len(ncol(X)), function(j) {
    others <- Rx[, -j, drop = FALSE]
    rx <- resid(lm(Rx[, j] ~ others)); ryr <- resid(lm(ry ~ others))
    suppressWarnings(cor(rx, ryr))
  })
}
res <- data.table(parameter = sens_pars)
for (out in colnames(Y)) res[[out]] <- round(prcc(X, Y[, out]), 3)
# approximate significance band for |PRCC| (n-dependent)
res$signif_band <- round(2 / sqrt(nrow(X)), 3)
setorder(res, -abs(soto_PFS))
fwrite(res, file.path(DIR_RESULTS, "sensitivity_prcc.csv"))
fwrite(as.data.table(cbind(X, Y)), file.path(DIR_RESULTS, "sensitivity_samples.csv"))
cat("\n== PRCC (|value| > signif_band is influential) ==\n"); print(res)
cat("Saved results/sensitivity_prcc.csv, sensitivity_samples.csv\n")
