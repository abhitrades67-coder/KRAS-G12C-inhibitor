#!/usr/bin/env Rscript
# 09_figures.R — all figures (PROTOCOL Phase 09). Reads the results/ artefacts produced by
# the earlier phases; each panel is guarded so the script runs with whatever is available.
# Usage: Rscript scripts/09_figures.R
source("R/setup.R"); source("R/calibrate.R")
suppressWarnings(suppressMessages({library(ggplot2); library(patchwork); library(scales); library(viridis)}))
theme_set(theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"), legend.position = "bottom"))
RES <- DIR_RESULTS; FIG <- DIR_FIGURES
sv  <- function(p, f, w = 8, h = 5.5) { ggsave(file.path(FIG, f), p, width = w, height = h, dpi = 150); cat("  figure:", f, "\n") }
arm_cols <- c(sotorasib = "#2c7fb8", adagrasib = "#7fcdbb", docetaxel = "#f03b20")
have <- function(f) file.exists(file.path(RES, f))
par <- load_fit()

# 1. Calibration: observed vs predicted ------------------------------------------------
if (have("calibration_targets_vs_pred.csv")) {
  d <- fread(file.path(RES, "calibration_targets_vs_pred.csv"))
  d[, grp := ifelse(metric == "occupancy", "Occupancy (fraction)", ifelse(metric=="ORR","ORR (%)","median PFS (months)"))]
  p <- ggplot(d, aes(value, pred, colour = arm)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = value - sd, xmax = value + sd), height = 0, alpha = .5) +
    geom_point(size = 3) + facet_wrap(~grp, scales = "free") +
    scale_colour_manual(values = arm_cols) +
    labs(title = "Calibration: model vs published anchors", x = "Observed", y = "Predicted", colour = NULL)
  sv(p, "fig01_calibration.png", 9, 4)
}

# 2. Waterfall (sotorasib best % change) ----------------------------------------------
if (have("vpop_outcomes.csv")) {
  o <- fread(file.path(RES, "vpop_outcomes.csv"))
  for (arm in intersect(ARMS, unique(o$arm))) {
    w <- o[arm == get("arm")][order(best_change)]; w[, x := .I]
    w[, resp := ifelse(best_change <= -30, "response", ifelse(best_change >= 20, "progression", "stable"))]
    p <- ggplot(w, aes(x, best_change, fill = resp)) + geom_col(width = 1) +
      geom_hline(yintercept = c(-30, 20), linetype = 2, colour = "grey40") +
      scale_fill_manual(values = c(response="#2c7fb8", stable="#bdbdbd", progression="#f03b20")) +
      coord_cartesian(ylim = c(-100, 100)) +
      labs(title = paste0("Best response waterfall — ", arm), x = "patients (sorted)",
           y = "best change in SLD (%)", fill = NULL)
    sv(p, sprintf("fig02_waterfall_%s.png", arm), 8, 4.5)
  }
}

# 3. Spider (SLD trajectories) --------------------------------------------------------
if (file.exists(file.path(RES, "vpop_trajectories.rds"))) {
  tr <- readRDS(file.path(RES, "vpop_trajectories.rds"))
  p <- ggplot(tr, aes(time/7, sld, group = interaction(arm, patient), colour = arm)) +
    geom_line(alpha = .25) + geom_hline(yintercept = 70, linetype = 2) +
    scale_colour_manual(values = arm_cols) + facet_wrap(~arm) +
    labs(title = "Tumour size trajectories (random subset)", x = "weeks", y = "SLD (% baseline)", colour = NULL) +
    coord_cartesian(ylim = c(0, 200))
  sv(p, "fig03_spider.png", 10, 4)
}

# 4. Kaplan-Meier PFS -----------------------------------------------------------------
if (have("vpop_outcomes.csv")) {
  o <- fread(file.path(RES, "vpop_outcomes.csv"))
  km <- rbindlist(lapply(unique(o$arm), function(a) {
    d <- o[arm == a]; s <- survfit(Surv(pfs_days, pfs_event) ~ 1, data = d)
    data.table(arm = a, time = c(0, s$time)/DAYS_PER_MONTH, surv = c(1, s$surv)) }))
  p <- ggplot(km, aes(time, 100*surv, colour = arm)) + geom_step(linewidth = 1) +
    scale_colour_manual(values = arm_cols) +
    labs(title = "Progression-free survival (in silico)", x = "months", y = "PFS (%)", colour = NULL) +
    coord_cartesian(xlim = c(0, 18))
  sv(p, "fig04_km_pfs.png", 8, 5)
}

# 5. Target engagement / signalling for a typical patient -----------------------------
{
  dt <- simulate_patient(par, arm_drive("sotorasib", par), fine = TRUE)
  m <- melt(dt[, .(week = time/7, `Covalent occupancy` = occ, `MAPK signal S` = S, `SLD/100` = SLD/100)],
            id.vars = "week")
  p <- ggplot(m, aes(week, value, colour = variable)) + geom_line(linewidth = 1) +
    scale_colour_manual(values = c("#2c7fb8","#d95f02","#1b9e77")) +
    labs(title = "Target engagement, pathway signal and tumour (typical sotorasib patient)",
         x = "weeks", y = "value", colour = NULL)
  sv(p, "fig05_target_engagement.png", 8, 5)
}

# 6. PRCC tornado ---------------------------------------------------------------------
if (have("sensitivity_prcc.csv")) {
  s <- fread(file.path(RES, "sensitivity_prcc.csv"))
  sm <- melt(s, id.vars = c("parameter","signif_band"), variable.name = "output", value.name = "prcc")
  p <- ggplot(sm, aes(reorder(parameter, abs(prcc)), prcc, fill = output)) +
    geom_col(position = "dodge") + coord_flip() +
    geom_hline(yintercept = c(-1,1)*s$signif_band[1], linetype = 2, colour = "grey50") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Global sensitivity (PRCC)", x = NULL, y = "partial rank correlation", fill = NULL)
  sv(p, "fig06_prcc.png", 8, 5)
}

# 7. Scenarios ------------------------------------------------------------------------
if (have("scenarios.csv")) {
  s <- fread(file.path(RES, "scenarios.csv")); s[, scenario := factor(scenario, levels = rev(scenario))]
  p1 <- ggplot(s, aes(scenario, ORR)) + geom_col(fill = "#2c7fb8") +
    geom_errorbar(aes(ymin = ORR_lo, ymax = ORR_hi), width = .3) + coord_flip() +
    labs(title = "ORR by scenario", x = NULL, y = "ORR (%)")
  p2 <- ggplot(s, aes(scenario, medPFS)) + geom_col(fill = "#1b9e77") + coord_flip() +
    labs(title = "Median PFS by scenario", x = NULL, y = "months")
  sv(p1 + p2, "fig07_scenarios.png", 13, 5)
}

# 8. In silico hazard ratios ----------------------------------------------------------
if (have("insilico_summary.csv")) {
  s <- fread(file.path(RES, "insilico_summary.csv"))
  p <- ggplot(s, aes(HR, comparison)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = HR_lo, xmax = HR_hi), height = .2) + geom_point(size = 3, colour = "#2c7fb8") +
    scale_x_log10() + labs(title = "In silico PFS hazard ratios (95% CrI)", x = "hazard ratio (log scale)", y = NULL)
  sv(p, "fig08_insilico_hr.png", 8, 3.5)
}

# 9. Benchmark (MBMA) -----------------------------------------------------------------
if (have("benchmark_onset.csv")) {
  b <- fread(file.path(RES, "benchmark_onset.csv"))
  p <- ggplot(melt(b, id.vars = "week", variable.name = "model", value.name = "frac"),
              aes(week, 100*frac, colour = model)) + geom_line(linewidth = 1) +
    scale_colour_manual(values = c(mbma = "#636363", mech = "#2c7fb8"),
                        labels = c("MBMA (Emax-in-time)","mechanistic QSP")) +
    labs(title = "Response-onset benchmark", x = "weeks", y = "cumulative responders (%)", colour = NULL)
  sv(p, "fig09_benchmark.png", 8, 4.5)
}

cat("Figures written to", FIG, "\n")
