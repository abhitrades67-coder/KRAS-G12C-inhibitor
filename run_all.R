# run_all.R — full KRAS G12C QSP pipeline. Each step runs in its own R process; output goes
# to logs/<step>.log. Run from the project root:
#   Rscript run_all.R              # all steps
#   Rscript run_all.R 05 06        # only steps whose names start with 05 or 06
# Set QSP_THREADS to control parallelism (default: all cores). Calibration is the long step.
ROOT <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), winslash = "/")
RSCRIPT <- file.path(R.home("bin"), "Rscript")
STEPS <- list(
  c("01_verify",        "scripts/01_verify.R"),
  c("02_calibrate",     "scripts/02_calibrate.R 160 700 4 160"),
  c("03_identifiability","scripts/03_identifiability.R"),
  c("05_vpop",          "scripts/05_vpop.R 1000"),
  c("05b_ensemble",     "scripts/05b_ensemble.R 1200 100 400 150"),
  c("06_insilico",      "scripts/06_insilico_trial.R 4000 300 500"),
  c("04_validate",      "scripts/04_validate.R"),
  c("07_sensitivity",   "scripts/07_sensitivity.R 220"),
  c("07b_structural",   "scripts/07b_structural.R 220"),
  c("55_benchmarks",    "scripts/55_benchmarks.R"),
  c("08_scenarios",     "scripts/08_scenarios.R 800"),
  c("09_figures",       "scripts/09_figures.R"),
  c("10_tables",        "scripts/10_tables.R")
)
sel <- commandArgs(trailingOnly = TRUE)
dir.create(file.path(ROOT, "logs"), showWarnings = FALSE)
for (s in STEPS) {
  if (length(sel) && !any(startsWith(s[1], sel))) next
  parts <- strsplit(s[2], " ")[[1]]
  log <- file.path(ROOT, "logs", paste0(s[1], ".log"))
  cat(sprintf("[%s] %-18s ... ", format(Sys.time(), "%H:%M:%S"), s[1]))
  t0 <- Sys.time()
  rc <- system2(RSCRIPT, c(shQuote(file.path(ROOT, parts[1])), parts[-1]), stdout = log, stderr = log)
  cat(sprintf("%s (%.1f min)\n", if (rc == 0) "ok" else paste("FAILED, see", log),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (rc != 0) stop("step failed: ", s[1])
}
cat("Pipeline complete. See results/, figures/, and REPORT.html (knit REPORT.Rmd).\n")
