# arms.R — map a trial arm (drug + regimen) to the model drive.
# pk_factor scales a patient's normalised exposure (1 = typical approved-dose exposure);
# dose_factor lets scenarios change the dose/schedule (affects exposure, hence engagement).

ARMS <- c("sotorasib", "adagrasib", "docetaxel")

arm_drive <- function(arm, par, pk_factor = 1, dose_factor = 1, kcov_factor = 1, doce_factor = 1) {
  cn <- pk_factor * dose_factor
  if (arm == "sotorasib") {
    list(cn = cn, kcov = par[["kcov_soto"]] * kcov_factor, edir = 0, kind = "g12c")
  } else if (arm == "adagrasib") {
    list(cn = cn, kcov = par[["kcov_adag"]] * kcov_factor, edir = 0, kind = "g12c")
  } else if (arm == "docetaxel") {
    # comparator: pulsatile cytotoxic kill, with per-patient chemo-sensitivity
    list(cn = 0, kcov = 0, edir = par[["kkill_doce"]] * doce_factor, kind = "cyto")
  } else stop("unknown arm: ", arm)
}

cat("arms.R loaded.\n")
