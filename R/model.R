# model.R — the QSP right-hand side and single-patient simulator.
#
# State (6 slow variables; sub-daily PK is averaged to steady-state exposure, which is
# appropriate for the monthly ORR/PFS endpoints and keeps the system non-stiff):
#   Rnc    non-covalent KRAS G12C protein      (baseline 1)
#   Rb     covalently drug-bound KRAS G12C     (baseline 0)
#   U      adaptive upstream input (RTK/SOS/SHP2, feedback-controlled, baseline 1)
#   res    acquired adaptive resistance in [0,1) (baseline 0) -> sensitive clone loses drug
#          suppression over weeks-to-months while engaged, driving acquired progression
#   Ts     sensitive tumour burden             (baseline 1-fr0)
#   Tr     resistant tumour burden             (baseline fr0)
#
# Drive per arm (a list `drv`):
#   cn       normalised free exposure of the G12C drug (1 = approved dose)
#   kcov     per-drug covalent capture rate at reference exposure (/day)
#   edir     direct cytotoxic (docetaxel) death rate (/day); full on the sensitive clone,
#            attenuated (x doce_res) on the drug-tolerant resistant clone
#
# Signal S = ppERK / ppERK_baseline (dimensionless) drives proliferation and,
# through its suppression, death. Feedback raises U when S < 1 (reactivation).

qsp_rhs <- function(t, y, p) {
  Rnc <- y[["Rnc"]]; Rb <- y[["Rb"]]; U <- y[["U"]]; res <- y[["res"]]
  Ts <- y[["Ts"]]; Tr <- y[["Tr"]]
  Rnc <- max(Rnc, 0); Ts <- max(Ts, 0); Tr <- max(Tr, 0); U <- max(U, 0); res <- clip(res, 0, 1)

  # active fraction among non-covalent protein (fast nucleotide-cycle QSS)
  a  <- (p$kact_base * U) / (p$kact_base * U + p$kgap)
  Rg <- (1 - a) * Rnc          # GDP-bound, drug-accessible
  Ra <- a * Rnc                # GTP-bound, active

  # signalling (normalised to baseline) ------------------------------------
  P   <- Ra / (Ra + p$Kd_sig)
  S   <- clip(P / p$P0, 0, 1.6)

  # adaptive feedback: upstream input relaxes toward target (up when ERK low)
  Utarget <- 1 + p$fb * max(0, 1 - S)
  dU <- p$kfb * (Utarget - U)

  # target engagement ------------------------------------------------------
  capture <- p$kcov_eff * p$cn * Rg            # irreversible covalent capture of GDP pool
  occ <- Rb / (Rnc + Rb + 1e-12)
  dRnc <- p$ksyn - p$kdeg * Rnc - capture
  dRb  <- capture - p$kdeg * Rb

  # acquired adaptive resistance (res in [0,1)): while the target is engaged, bypass/amplification
  # accumulates and the sensitive clone's EFFECTIVE pathway signal returns toward baseline, so the
  # drug's hold wanes over weeks-to-months and the tumour re-grows. kprog sets the timescale (PFS).
  dres <- p$kprog * occ * (1 - res)

  # direct cytotoxic (docetaxel) effect (0 for G12C arms), waning over time as chemoresistance
  # develops (short, non-durable benefit -> shorter PFS than targeted therapy).
  edir <- p$edir_amp * exp(-p$kprog_doce * t)

  # tumour: two clones. The sensitive clone is KRAS-addicted (dependence dep) and ACQUIRES
  # resistance over time (its effective signal Ss returns to baseline as res->1). The resistant
  # clone carries PRE-EXISTING bypass (effective signal Sr stays high) and is only weakly touched
  # by chemotherapy. Growth = kprol*(d*Se + 1-d) - [kdeath0 + kkill*d*(1-Se) + kill].
  d  <- p$dep
  Ss <- clip(S + (1 - S) * res,      0, 1.6)   # sensitive: suppression lost as res grows
  Sr <- clip(S + (1 - S) * p$bypass, 0, 1.6)   # resistant: pre-existing bypass
  gS <- p$kprol * (d * Ss + (1 - d)) - (p$kdeath0 + p$kkill * d * (1 - Ss) + edir)
  gR <- p$kprol * (d * Sr + (1 - d)) - (p$kdeath0 + p$kkill * d * (1 - Sr) + edir * p$doce_res)
  dTs <- gS * Ts - p$kres * Ts
  dTr <- gR * Tr + p$kres * Ts

  list(c(Rnc = dRnc, Rb = dRb, U = dU, res = dres, Ts = dTs, Tr = dTr),
       S = S, Ss = Ss, a = a, edir = edir, occ = occ, res = res)
}

# Build a fully-resolved parameter environment for one patient/arm.
# `par` is a named numeric vector (PAR0 or a perturbed copy); `drv` sets the arm drive.
make_p <- function(par, drv) {
  p <- as.list(par)
  p$ksyn <- p$kdeg                     # baseline total protein = 1
  a0     <- (p$kact_base * 1) / (p$kact_base * 1 + p$kgap)
  Ra0    <- a0 * 1
  p$P0   <- Ra0 / (Ra0 + p$Kd_sig)     # baseline signal normaliser
  p$kcov_eff <- drv$kcov %||% 0
  p$cn       <- drv$cn   %||% 0
  p$edir_amp <- drv$edir %||% 0          # mean daily docetaxel kill (0 for G12C arms)
  p
}

# initial state at baseline (pre-treatment steady state)
init_state <- function(par) {
  fr0 <- clip(par[["fr0"]], 1e-4, 0.95)
  c(Rnc = 1, Rb = 0, U = 1, res = 0, Ts = 1 - fr0, Tr = fr0)
}

# Simulate one patient on one arm. Returns a data.table at the requested times (default the
# assessment grid; fine=TRUE returns the dense internal grid) with SLD (% baseline), the
# pathway signal S, the active fraction a, and covalent occupancy occ.
simulate_patient <- function(par, drv, times = ASSESS_GRID, fine = FALSE) {
  p  <- make_p(par, drv)
  y0 <- init_state(par)
  # all arms are smooth: coarse grid, denser early to resolve initial response kinetics.
  tt <- sort(unique(c(seq(0, min(126, max(times)), by = 7), seq(0, max(times), by = 21), times)))
  out <- tryCatch(
    as.data.frame(ode(y0, tt, function(t, y, pp) qsp_rhs(t, y, pp), p,
                      method = "lsoda", rtol = 1e-6, atol = 1e-8)),
    error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out$Ts)) || any(!is.finite(out$Tr))) return(NULL)
  dt <- as.data.table(out)
  dt[, SLD := 100 * (Ts + Tr)]                 # % of baseline (baseline Ts+Tr = 1)
  if (!fine) dt <- dt[time %in% times]
  dt[]
}

cat("model.R loaded.\n")
