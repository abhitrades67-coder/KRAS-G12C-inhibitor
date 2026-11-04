#!/usr/bin/env Rscript
# 03_identifiability.R — practical identifiability of the fitted parameters (PROTOCOL Phase 03):
#   * Hessian-based standard errors and the parameter correlation matrix (confounded pairs)
#   * prior->posterior shrinkage (posterior 95% width vs prior bound width)
#   * a conditional profile of the objective around the optimum for each parameter
# Usage: Rscript scripts/03_identifiability.R
source("R/setup.R"); source("R/calibrate.R")
suppressWarnings(suppressMessages(library(numDeriv)))
fit <- readRDS(file.path(DIR_RESULTS, "fit.rds")); best <- fit$best
N <- 150
cat("Identifiability at the fitted optimum (N =", N, ")\n")

# objective as negative log-likelihood: g = sum w z^2 = 2*NLL  =>  cov ~ 2 * inv(Hessian(g))
g <- function(x) { names(x) <- FIT_PARS; calib_objective(x, N = N, cores = QSP_THREADS) }
cat("Computing Hessian (may take a few minutes)...\n")
H <- tryCatch(hessian(g, best, method.args = list(eps = 1e-3, d = 0.02, r = 2)),
              error = function(e) NULL)

se <- rep(NA_real_, length(best)); cormat <- NULL
if (!is.null(H) && all(is.finite(H))) {
  Hs <- (H + t(H)) / 2
  cov <- tryCatch(2 * solve(Hs), error = function(e) 2 * MASS::ginv(Hs))
  d <- diag(cov); se <- ifelse(d > 0, sqrt(d), NA_real_)
  sdinv <- ifelse(is.na(se), 0, 1/se); cormat <- diag(sdinv) %*% cov %*% diag(sdinv)
  dimnames(cormat) <- list(FIT_PARS, FIT_PARS)
}

prior_w <- (FIT_UB - FIT_LB)[FIT_PARS]
post_w  <- 2 * 1.96 * se                       # approx 95% posterior width
shrink  <- 1 - post_w / prior_w
idtab <- data.table(parameter = FIT_PARS, estimate = as.numeric(best),
                    se = se, post95_width = post_w, prior_width = as.numeric(prior_w),
                    shrinkage = shrink,
                    identified = ifelse(is.na(se), "no",
                                 ifelse(shrink > 0.5, "yes", ifelse(shrink > 0.2, "partial", "weak"))))
setorder(idtab, -shrinkage)
fwrite(idtab, file.path(DIR_RESULTS, "identifiability.csv"))
cat("\n== Prior->posterior shrinkage ==\n"); print(idtab[, .(parameter, estimate = round(estimate,4),
    se = round(se,4), shrinkage = round(shrinkage,2), identified)])

if (!is.null(cormat)) {
  fwrite(as.data.table(round(cormat, 3), keep.rownames = "param"),
         file.path(DIR_RESULTS, "identifiability_correlation.csv"))
  pr <- which(abs(cormat) > 0.7 & upper.tri(cormat), arr.ind = TRUE)
  cat("\n== Confounded parameter pairs (|corr| > 0.7) ==\n")
  if (nrow(pr)) for (k in seq_len(nrow(pr)))
    cat(sprintf("  %s <-> %s : %.2f\n", FIT_PARS[pr[k,1]], FIT_PARS[pr[k,2]], cormat[pr[k,1], pr[k,2]]))
  else cat("  none\n")
}

# conditional profile (fix one parameter, hold others at optimum) around +/-40%
cat("\n== Conditional objective profiles ==\n")
prof <- list()
for (p in FIT_PARS) {
  grid <- seq(0.6, 1.4, by = 0.1) * best[p]
  grid <- pmin(FIT_UB[p], pmax(FIT_LB[p], grid))
  gv <- vapply(grid, function(v) { x <- best; x[p] <- v; g(x) }, numeric(1))
  prof[[p]] <- data.table(parameter = p, value = grid, objective = gv)
}
profile_dt <- rbindlist(prof)
fwrite(profile_dt, file.path(DIR_RESULTS, "identifiability_profiles.csv"))

# structural-identifiability notes (a priori, honest statement of what is/ isn't separable)
notes <- c(
 "Structural / a-priori identifiability notes:",
 "- Covalent potency (kcov_*) and free-exposure scale are only jointly identified; the exposure",
 "  scale (Cref) is fixed from PK so kcov_* carries the engagement efficiency.",
 "- Effect-site partition is not identifiable from these endpoints; effect-site = free plasma (fixed).",
 "- kprol and kdeath0 trade off in the net growth rate; only the net rate is informed, so both are",
 "  fixed from natural-history doubling time and only drug/resistance parameters are estimated.",
 "- A full differential-algebra structural analysis is out of scope; practical identifiability is",
 "  assessed here by the Hessian covariance, prior->posterior shrinkage and the objective profiles.")
writeLines(notes, file.path(DIR_RESULTS, "identifiability_notes.txt"))
cat(paste(notes, collapse = "\n"), "\n")
cat("\nSaved results/identifiability*.csv, identifiability_notes.txt\n")
