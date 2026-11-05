#!/usr/bin/env Rscript
# 07b_structural.R — structural model qualification (PROTOCOL Phase 7, structure).
# Refit reduced variants of the model (each disabling one mechanism) and compare their fit to
# the full model by objective and an AIC-like score (objective + 2*n_free). If removing a
# mechanism clearly worsens the fit, that mechanism is warranted. Usage: Rscript scripts/07b_structural.R [N]
source("R/setup.R"); source("R/calibrate.R")
a <- as.numeric(commandArgs(trailingOnly = TRUE)); N <- if (length(a) && !is.na(a[1])) a[1] else 220
fit <- readRDS(file.path(DIR_RESULTS, "fit.rds")); best <- fit$best

# each variant fixes some parameters (disabling a mechanism) and refits the rest
variants <- list(
  list(name = "full model",                 fix = c()),
  list(name = "no adaptive feedback (fb=0)", fix = c(fb = 0)),
  list(name = "no acquired resistance (kprog=0)", fix = c(kprog = 0)),
  list(name = "no pre-existing resistance (fr0~0)", fix = c(fr0 = 0.02)),
  list(name = "no feedback & no acquired res", fix = c(fb = 0, kprog = 0))
)

refit_fix <- function(fixv) {
  free <- setdiff(FIT_PARS, names(fixv))
  th0  <- best; if (length(fixv)) th0[names(fixv)] <- fixv
  fn <- function(x) { th <- th0; th[free] <- pmin(FIT_UB[free], pmax(FIT_LB[free], x))
                      calib_objective(th, N = N, seed = SEED, cores = QSP_THREADS) }
  op <- optim(th0[free], fn, method = "Nelder-Mead", control = list(maxit = 120, reltol = 1e-4))
  th <- th0; th[free] <- pmin(FIT_UB[free], pmax(FIT_LB[free], op$par))
  list(theta = th, obj = op$value, k = length(free))
}

rows <- list()
for (v in variants) {
  r <- refit_fix(v$fix)
  pred <- predict_endpoints(assemble_par(r$theta), N = 400, seed = SEED + 1, cores = QSP_THREADS)
  obj_eval <- sum(vapply(seq_len(nrow(TARGETS)), function(i)
    TARGETS$w[i] * ((pred[[TARGETS$tkey[i]]] - TARGETS$value[i]) / TARGETS$sd[i])^2, numeric(1)))
  rows[[v$name]] <- data.table(model = v$name, n_free = r$k, objective = round(obj_eval, 2),
    AIC_like = round(obj_eval + 2 * r$k, 2),
    soto_ORR = round(pred[["soto_ORR"]],1), soto_PFS = round(pred[["soto_PFS"]],1),
    adag_ORR = round(pred[["adag_ORR"]],1), adag_PFS = round(pred[["adag_PFS"]],1),
    doce_ORR = round(pred[["doce_ORR"]],1), doce_PFS = round(pred[["doce_PFS"]],1))
  cat(sprintf("  %-34s obj=%.2f AIC~%.2f\n", v$name, obj_eval, obj_eval + 2*r$k))
}
st <- rbindlist(rows); setorder(st, AIC_like)
st[, dAIC := round(AIC_like - min(AIC_like), 2)]
fwrite(st, file.path(DIR_RESULTS, "structural_comparison.csv"))
cat("\n== Structural comparison (lower AIC-like = better; dAIC vs best) ==\n")
print(st[, .(model, n_free, objective, AIC_like, dAIC)])
cat("Saved results/structural_comparison.csv\n")
