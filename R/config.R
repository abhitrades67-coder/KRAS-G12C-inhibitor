# config.R — paths, global constants, seeds, small helpers.
# Everything downstream sources this first (directly or via 00_setup.R).

suppressWarnings(suppressMessages({
  library(deSolve)
  library(data.table)
}))

# --- project root (works whether sourced from root or from scripts/) ---
PROJ <- tryCatch({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(file.path(dirname(sub("^--file=", "", a[1])), ".."), winslash = "/", mustWork = FALSE)
  else getwd()
}, error = function(e) getwd())
# If sourced interactively from the project root, prefer the working directory when it has R/.
if (!dir.exists(file.path(PROJ, "R")) && dir.exists(file.path(getwd(), "R"))) PROJ <- getwd()
if (basename(PROJ) == "scripts") PROJ <- dirname(PROJ)

DIR_R        <- file.path(PROJ, "R")
DIR_SCRIPTS  <- file.path(PROJ, "scripts")
DIR_RESULTS  <- file.path(PROJ, "results")
DIR_FIGURES  <- file.path(PROJ, "figures")
DIR_LOGS     <- file.path(PROJ, "logs")
DIR_CALSTART <- file.path(PROJ, "calibration_start")
for (d in c(DIR_RESULTS, DIR_FIGURES, DIR_LOGS, DIR_CALSTART)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# --- reproducibility & compute ---
SEED <- 20261104L
set.seed(SEED)
QSP_THREADS <- as.integer(Sys.getenv("QSP_THREADS", unset = NA))
if (is.na(QSP_THREADS)) QSP_THREADS <- max(1L, parallel::detectCores())

# --- simulation clock (days) ---
HORIZON_D     <- 540           # ~18 months follow-up
ASSESS_EVERY  <- 42            # RECIST assessment interval (6 weeks), CodeBreaK/KRYSTAL-like
ASSESS_GRID   <- seq(0, HORIZON_D, by = ASSESS_EVERY)

# --- RECIST thresholds (on sum of longest diameters, SLD, % of baseline) ---
RECIST_PR <- 70    # <=70% of baseline  => >=30% reduction  => (partial) response
RECIST_PD_REL <- 1.20  # >=20% increase over nadir
RECIST_PD_ABS <- 5     # plus a small absolute increase (percentage points of baseline)

cat(sprintf("config.R loaded. PROJ=%s  threads=%d  seed=%d\n", PROJ, QSP_THREADS, SEED))

# small helpers ---------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a)==1 && is.na(a))) b else a
clip   <- function(x, lo, hi) pmin(hi, pmax(lo, x))
logit  <- function(p) log(p/(1-p)); invlogit <- function(x) 1/(1+exp(-x))
