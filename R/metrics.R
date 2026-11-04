# metrics.R — RECIST response & PFS from an SLD time-course, virtual-population
# generation, and per-arm population summaries (ORR %, median PFS, occupancy).

suppressWarnings(suppressMessages(library(survival)))
DAYS_PER_MONTH <- 30.4375

# confirmed objective response (1 = PR/CR): SLD <= RECIST_PR % of baseline at two
# consecutive assessments (RECIST confirmation). This keeps durable responses and
# discards the transient between-cycle dips of pulsatile chemotherapy.
best_response <- function(sld) {
  resp <- sld <= RECIST_PR
  if (length(resp) >= 2 && any(resp[-length(resp)] & resp[-1], na.rm = TRUE)) return(1L)
  # allow an unconfirmable response only if it is the final assessment and sustained
  0L
}

# time to RECIST progression on the assessment grid (days); NA if none by horizon
pfs_event <- function(times, sld) {
  nadir <- sld[1]                       # starts at baseline (100)
  for (i in seq_along(times)[-1]) {
    nadir <- min(nadir, sld[i - 1])
    if (sld[i] >= RECIST_PD_REL * nadir && sld[i] >= nadir + RECIST_PD_ABS)
      return(times[i])
  }
  NA_real_
}

# ---- virtual population: deterministic given seed (common random numbers) ----
gen_vpop <- function(par, N, seed = SEED) {
  vars <- names(POP_CV)
  set.seed(seed)
  Z <- matrix(rnorm(N * length(vars)), N, length(vars), dimnames = list(NULL, vars))
  lognf <- function(cv, z) exp(cv * z - 0.5 * cv^2)   # median-preserving lognormal factor
  pkZ <- rnorm(N)                                     # independent exposure draw
  dcZ <- rnorm(N)                                     # independent docetaxel-sensitivity draw
  lapply(seq_len(N), function(i) {
    pi <- par
    # KRAS-pathway dependence on the logit scale -> proper (0,1) spread (main ORR driver)
    pi["dep"]    <- invlogit(logit(clip(par["dep"],0.02,0.98)) + POP_CV["dep"] * Z[i,"dep"])
    pi["kkill"]  <- par["kkill"]  * lognf(POP_CV["kkill"],  Z[i,"kkill"])
    pi["bypass"] <- clip(par["bypass"] * lognf(POP_CV["bypass"], Z[i,"bypass"]), 0.01, 0.99)
    pi["fr0"]    <- clip(par["fr0"]    * lognf(POP_CV["fr0"],    Z[i,"fr0"]),    1e-4, 0.9)
    pi["kprol"]  <- par["kprol"]  * lognf(POP_CV["kprol"],  Z[i,"kprol"])
    pi["fb"]     <- par["fb"]     * lognf(POP_CV["fb"],     Z[i,"fb"])
    list(par = pi,
         pk_factor   = lognf(0.25, pkZ[i]),             # G12C exposure heterogeneity
         kcov_factor = lognf(POP_CV["kcov"], Z[i,"kcov"]),
         doce_factor = lognf(0.20, dcZ[i]))             # docetaxel chemo-sensitivity heterogeneity
  })
}

# simulate one arm across a virtual population; returns per-patient SLD + summary
simulate_arm <- function(par, arm, vpop, cores = QSP_THREADS) {
  f <- function(v) {
    drv <- arm_drive(arm, v$par, pk_factor = v$pk_factor, kcov_factor = v$kcov_factor,
                     doce_factor = v$doce_factor %||% 1)
    dt  <- simulate_patient(v$par, drv)
    if (is.null(dt)) return(NULL)
    list(sld = dt$SLD, times = dt$time, occ = dt$occ, S = dt$S)
  }
  res <- if (cores > 1) parallel::mclapply(vpop, f, mc.cores = cores) else lapply(vpop, f)
  ok  <- !vapply(res, is.null, logical(1))
  res[ok]
}

# population summary for one arm
summarise_arm <- function(par, arm, vpop, cores = QSP_THREADS) {
  res <- simulate_arm(par, arm, vpop, cores)
  n <- length(res); if (n == 0) return(list(ORR = NA, medPFS = NA, occ = NA, n = 0))
  resp <- vapply(res, function(r) best_response(r$sld), numeric(1))
  pfs  <- vapply(res, function(r) pfs_event(r$times, r$sld), numeric(1))
  ev   <- as.integer(!is.na(pfs))
  tt   <- ifelse(is.na(pfs), HORIZON_D, pfs)
  med  <- tryCatch({
    km <- survfit(Surv(tt, ev) ~ 1)
    unname(summary(km)$table["median"]) / DAYS_PER_MONTH
  }, error = function(e) NA_real_)
  if (is.na(med)) med <- median(tt) / DAYS_PER_MONTH    # fallback if KM median not reached
  occ_plateau <- mean(vapply(res, function(r) {
    j <- which.min(abs(r$times - 84)); r$occ[j]          # ~12-week plateau occupancy
  }, numeric(1)))
  list(ORR = 100 * mean(resp), medPFS = med, occ = occ_plateau, n = n,
       pfs_days = tt, pfs_event = ev, resp = resp)
}

cat("metrics.R loaded.\n")
