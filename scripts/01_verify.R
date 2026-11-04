#!/usr/bin/env Rscript
# 01_verify.R — numerical & structural verification of the model (PROTOCOL Phase 01).
# Fast, deterministic checks; writes results/verification.csv and stops on hard failures.
source("R/setup.R")
chk <- list(); add <- function(name, pass, detail = "") {
  chk[[length(chk) + 1]] <<- data.table(check = name, pass = pass, detail = detail)
  cat(sprintf("[%s] %s %s\n", if (pass) "PASS" else "FAIL", name, detail))
}

par <- PAR0

# 1. Baseline steady state: with no drug, protein total = 1, occupancy = 0, signal S = 1.
drv0 <- list(cn = 0, kcov = 0, edir = 0)
d0 <- simulate_patient(par, drv0, times = c(0, 90, 180), fine = FALSE)
add("baseline protein ~1",  abs((d0$Rnc[1] + d0$Rb[1]) - 1) < 1e-6, sprintf("=%.4f", d0$Rnc[1] + d0$Rb[1]))
add("baseline occupancy 0", d0$occ[1] < 1e-8, sprintf("=%.2e", d0$occ[1]))
add("baseline signal S~1",  abs(d0$S[1] - 1) < 1e-6, sprintf("=%.5f", d0$S[1]))

# 2. Untreated tumour grows at the intended net rate (kprol-kdeath0) and is monotone.
net <- par[["kprol"]] - par[["kdeath0"]]
growth_ok <- all(diff(d0$SLD) > 0) && abs(log(d0$SLD[3]/d0$SLD[1]) / 180 - net) < 0.002
add("untreated growth = net rate", growth_ok, sprintf("net=%.4f/day", net))

# 3. Mass/positivity: all states finite and non-negative over a full treated run.
dt <- simulate_patient(par, arm_drive("sotorasib", par), fine = TRUE)
add("states finite & non-negative",
    all(is.finite(as.matrix(dt[, .(Rnc, Rb, U, Ts, Tr)]))) &&
    min(dt$Rnc, dt$Rb, dt$Ts, dt$Tr) >= -1e-8, "")

# 4. Monotonic occupancy in exposure: higher kcov -> higher plateau occupancy.
occ_of <- function(k) { p <- par; p["kcov_soto"] <- k
  tail(simulate_patient(p, arm_drive("sotorasib", p), times = c(0, 120))$occ, 1) }
o_lo <- occ_of(3); o_hi <- occ_of(20)
add("occupancy increases with exposure", o_hi > o_lo, sprintf("%.2f -> %.2f", o_lo, o_hi))

# 5. Dose-response monotonicity: deeper signal suppression -> deeper nadir (lower best SLD).
nadir_of <- function(dose) { d <- arm_drive("sotorasib", par, dose_factor = dose)
  min(simulate_patient(par, d, fine = TRUE)$SLD) }
add("nadir deepens with dose", nadir_of(2) <= nadir_of(0.3) + 1e-6,
    sprintf("%.1f vs %.1f", nadir_of(2), nadir_of(0.3)))

# 6. Feedback raises U under drug (adaptive reactivation present).
add("feedback raises U under drug", tail(dt$U, 1) > 1 + 1e-3, sprintf("U=%.2f", tail(dt$U, 1)))

# 7. Integrator consistency: tighter tolerance changes endpoint SLD negligibly.
p <- par
loose <- simulate_patient(p, arm_drive("sotorasib", p), times = c(0, 180))$SLD[2]
f_tight <- function() { environment(qsp_rhs); as.data.frame(deSolve::ode(
  init_state(p), c(0,180), function(t,y,pp) qsp_rhs(t,y,pp), make_p(p, arm_drive("sotorasib", p)),
  method = "lsoda", rtol = 1e-9, atol = 1e-11)) }
tight <- 100 * (f_tight()$Ts[2] + f_tight()$Tr[2])
add("integrator tolerance-stable", abs(loose - tight) < 0.5, sprintf("|d|=%.3f", abs(loose - tight)))

# 8. PFS/response logic unit tests on synthetic series.
add("response logic: deep+confirmed = PR",
    best_response(c(100, 50, 45, 40)) == 1 && best_response(c(100, 60, 90, 95)) == 0, "")
add("PFS logic: +20% over nadir triggers",
    pfs_event(c(0,42,84,126), c(100, 50, 55, 62)) == 126 &&
    is.na(pfs_event(c(0,42,84), c(100, 80, 70))), "")

res <- rbindlist(chk)
fwrite(res, file.path(DIR_RESULTS, "verification.csv"))
cat(sprintf("\n%d/%d checks passed. Saved results/verification.csv\n", sum(res$pass), nrow(res)))
if (!all(res$pass)) quit(status = 1)
