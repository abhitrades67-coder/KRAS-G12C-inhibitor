# calibrate.R — map the fitted-parameter vector to predicted endpoints and an objective.
# Used by scripts/02_calibrate.R and by the identifiability / validation scripts.

# Build a full parameter vector from a fitted-subset vector theta (named by FIT_PARS).
assemble_par <- function(theta, base = PAR0) {
  p <- base
  p[names(theta)] <- as.numeric(theta)
  p
}

# Predict the model endpoints for a parameter set, as a named vector keyed like TARGETS$tkey.
predict_endpoints <- function(par, N = 160, seed = SEED, cores = 1) {
  vp   <- gen_vpop(par, N, seed)
  soto <- summarise_arm(par, "sotorasib", vp, cores)
  adag <- summarise_arm(par, "adagrasib", vp, cores)
  doce <- summarise_arm(par, "docetaxel", vp, cores)
  c(soto_ORR = soto$ORR, soto_PFS = soto$medPFS, soto_occ = soto$occ,
    adag_ORR = adag$ORR, adag_PFS = adag$medPFS,
    doce_ORR = doce$ORR, doce_PFS = doce$medPFS)
}

# Weighted sum-of-squared standardised residuals against TARGETS.
calib_objective <- function(theta_named, N = 160, seed = SEED, cores = 1, verbose = FALSE) {
  theta_named[] <- pmin(FIT_UB[names(theta_named)], pmax(FIT_LB[names(theta_named)], theta_named))
  par  <- assemble_par(theta_named)
  pred <- tryCatch(predict_endpoints(par, N, seed, cores), error = function(e) NULL)
  if (is.null(pred)) return(1e6)
  obj <- 0
  for (i in seq_len(nrow(TARGETS))) {
    k <- TARGETS$tkey[i]; pv <- pred[[k]]
    if (is.null(pv) || !is.finite(pv)) { obj <- obj + 50; next }
    obj <- obj + TARGETS$w[i] * ((pv - TARGETS$value[i]) / TARGETS$sd[i])^2
  }
  if (verbose) cat(sprintf("obj=%.2f | %s\n", obj,
      paste(sprintf("%s=%.1f", TARGETS$tkey, pred[TARGETS$tkey]), collapse=" ")))
  obj
}

cat("calibrate.R loaded.\n")

# Load the fitted parameter vector (falls back to PAR0 if calibration has not run yet).
load_fit <- function() {
  fp <- file.path(DIR_RESULTS, "fit.rds")
  if (file.exists(fp)) readRDS(fp)$par_fit else { message("fit.rds not found; using PAR0"); PAR0 }
}

# Objective restricted to a subset of target keys (used for leave-one-drug-out validation).
calib_objective_sub <- function(theta_named, keys, N = 140, seed = SEED, cores = QSP_THREADS) {
  theta_named[] <- pmin(FIT_UB[names(theta_named)], pmax(FIT_LB[names(theta_named)], theta_named))
  par  <- assemble_par(theta_named)
  pred <- tryCatch(predict_endpoints(par, N, seed, cores), error = function(e) NULL)
  if (is.null(pred)) return(1e6)
  tt <- TARGETS[tkey %in% keys]
  sum(vapply(seq_len(nrow(tt)), function(i) {
    pv <- pred[[tt$tkey[i]]]; if (is.null(pv) || !is.finite(pv)) return(50)
    tt$w[i] * ((pv - tt$value[i]) / tt$sd[i])^2 }, numeric(1)))
}
