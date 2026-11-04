#!/usr/bin/env Rscript
# 04_validate.R — out-of-sample validation (PROTOCOL Phase 06 strategy):
#   LODO  leave-one-drug-out: refit the SHARED pathway/resistance parameters on the other
#         drugs (each drug's own covalent/cytotoxic potency is treated as preclinically
#         transferred and held fixed), then predict the withheld drug's ORR and PFS.
#   TEMPORAL: fit to the response (early) endpoints only and predict PFS (late) endpoint.
# Pre-specified acceptance: |ORR| within 10 pp and |PFS| within 1.5 months.
# Usage: Rscript scripts/04_validate.R
source("R/setup.R"); source("R/calibrate.R")
fit <- readRDS(file.path(DIR_RESULTS, "fit.rds")); best <- fit$best
N <- 180; Nval <- 400
shared <- c("fb","kkill","dep","bypass","fr0","kres")
ownpar <- list(sotorasib = "kcov_soto", adagrasib = "kcov_adag", docetaxel = "kkill_doce")
dkeys  <- list(sotorasib = c("soto_ORR","soto_PFS","soto_occ"),
               adagrasib = c("adag_ORR","adag_PFS"),
               docetaxel = c("doce_ORR","doce_PFS"))

refit <- function(free, keys, start) {
  fn <- function(x) { th <- best; th[free] <- pmin(FIT_UB[free], pmax(FIT_LB[free], x))
                      calib_objective_sub(th, keys, N = N) }
  op <- optim(start[free], fn, method = "Nelder-Mead", control = list(maxit = 80, reltol = 1e-4))
  th <- best; th[free] <- pmin(FIT_UB[free], pmax(FIT_LB[free], op$par)); th
}

cat("== LODO (leave-one-drug-out): predict a withheld drug from the shared model ==\n")
lodo <- list()
for (d in names(dkeys)) {
  keep_keys <- TARGETS$tkey[!TARGETS$tkey %in% dkeys[[d]]]
  free <- c(shared, unlist(ownpar[setdiff(names(ownpar), d)]))   # drop the withheld drug's own potency
  th   <- refit(free, keep_keys, best)                          # withheld drug's own potency stays = best
  par  <- assemble_par(th)
  pred <- predict_endpoints(par, N = Nval, seed = SEED + 7, cores = QSP_THREADS)
  for (k in dkeys[[d]]) {
    obs <- TARGETS[tkey == k]$value
    lodo[[length(lodo)+1]] <- data.table(scheme = "LODO", held_out = d, key = k,
      metric = TARGETS[tkey==k]$metric, observed = obs, predicted = round(pred[[k]],2),
      abs_err = round(abs(pred[[k]] - obs), 2))
  }
  cat(sprintf("  withheld %-10s -> %s\n", d,
      paste(sprintf("%s obs=%.1f pred=%.1f", dkeys[[d]], TARGETS[tkey %in% dkeys[[d]]]$value,
                    pred[dkeys[[d]]]), collapse = " | ")))
}
lodo_dt <- rbindlist(lodo)

cat("\n== TEMPORAL: fit response (ORR+occupancy) only, predict PFS ==\n")
early <- TARGETS$tkey[TARGETS$metric %in% c("ORR","occupancy")]
th_t  <- refit(FIT_PARS, early, best)
par_t <- assemble_par(th_t)
pred_t<- predict_endpoints(par_t, N = Nval, seed = SEED + 8, cores = QSP_THREADS)
temp <- rbindlist(lapply(TARGETS[metric == "PFS"]$tkey, function(k) {
  obs <- TARGETS[tkey==k]$value
  data.table(scheme = "TEMPORAL(fit ORR->predict PFS)", held_out = "all PFS", key = k,
             metric = "PFS", observed = obs, predicted = round(pred_t[[k]],2),
             abs_err = round(abs(pred_t[[k]] - obs),2)) }))
print(temp[, .(key, observed, predicted, abs_err)])

val <- rbind(lodo_dt, temp)
val[, pass := ifelse(metric == "ORR", abs_err <= 10, ifelse(metric == "PFS", abs_err <= 1.5, abs_err <= 0.1))]
fwrite(val, file.path(DIR_RESULTS, "validation.csv"))
cat(sprintf("\nValidation: %d/%d held-out predictions within the pre-specified tolerance.\n",
            sum(val$pass, na.rm = TRUE), nrow(val)))
cat("Saved results/validation.csv\n")
