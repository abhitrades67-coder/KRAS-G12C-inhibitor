# params.R — parameters, per-drug PK, priors, and calibration targets.
# All rates are per DAY. Protein pools are normalised to baseline total = 1.
# Sources for every clinically anchored number are tagged here and expanded in
# verification/SOURCE_VERIFICATION.md.

# ---------------------------------------------------------------------------
# 1. Per-drug pharmacokinetics (one-compartment oral, apparent V/F and CL/F).
#    Concentrations are normalised to each drug's own reference exposure (Cref),
#    so the covalent-capture rate constant is a per-day rate at typical exposure.
# ---------------------------------------------------------------------------
PK <- list(
  sotorasib = list(                 # AMG 510, 960 mg once daily
    mw   = 560.6,                    # g/mol
    dose_mg = 960, interval_d = 1,   # QD
    ka   = 24.0,                     # /day  (rapid oral absorption, Tmax ~1 h)
    V    = 230,                      # L, apparent Vz/F (label range ~210-280 L)
    CL   = 32,                       # L/day? -> set via t1/2; see below
    thalf_h = 5.0,                   # terminal half-life ~5 h (sotorasib label)
    fu   = 0.11,                     # fraction unbound (~89% protein bound)
    kind = "g12c"
  ),
  adagrasib = list(                  # MRTX849, 600 mg twice daily
    mw   = 604.1,
    dose_mg = 600, interval_d = 0.5, # BID
    ka   = 12.0,                     # /day (Tmax ~6 h)
    V    = 1100,                     # L, large apparent volume
    thalf_h = 23.0,                  # long half-life ~23 h
    fu   = 0.02,                     # highly protein bound
    kind = "g12c"
  ),
  docetaxel = list(                  # 75 mg/m^2 IV every 3 weeks (comparator)
    mw   = 807.9,
    dose_mg = 140, interval_d = 21,  # ~75 mg/m^2 * 1.86 m^2 ~ 140 mg
    ka   = NA,                       # IV bolus into central
    V    = 75,                       # L (docetaxel Vss ~ 80-110 L; apparent)
    thalf_h = 11.0,                  # effective beta half-life (hours)
    fu   = 0.06,
    kind = "cyto"
  )
)
# derive elimination rate kel (/day) and a reference free concentration Cref (nM)
for (nm in names(PK)) {
  p <- PK[[nm]]
  kel <- log(2) / (p$thalf_h/24)            # /day
  PK[[nm]]$kel <- kel
  dose_nmol <- p$dose_mg / p$mw * 1e6        # nmol
  # average total concentration at steady state ~ (F*Dose)/(CL*tau); use kel*V as CL
  Css_avg_total <- dose_nmol / (kel * p$V * p$interval_d)   # nM (apparent)
  PK[[nm]]$Cref <- max(1e-6, Css_avg_total * p$fu)          # free reference, nM
}

# ---------------------------------------------------------------------------
# 2. Mechanistic parameters (population means). Fitted ones flagged in FIT_PARS.
# ---------------------------------------------------------------------------
PAR0 <- c(
  # --- target & cycling (mostly fixed from structure) ---
  kdeg      = 0.7,      # /day, KRAS G12C protein turnover (t1/2 ~1 day)
  kact_base = 16.0,     # /day, baseline nucleotide loading rate scale
  kgap      = 24.0,     # /day, intrinsic + GAP-stimulated hydrolysis
  # => baseline active fraction a0 = 16/(16+24) = 0.40
  Kd_sig    = 0.25,     # signalling half-saturation (in active-pool units)
  # --- adaptive feedback ---
  fb        = 2.6,      # FIT: feedback strength (ERK-relief -> upstream input)
  kfb       = 0.6,      # /day, feedback adaptation rate (reactivation over days)
  kprog     = 0.007,    # FIT: /day acquired adaptive-resistance rate (sets duration of benefit)
  # --- covalent capture efficiency (per drug), at reference exposure ---
  kcov_soto = 14.0,      # FIT: /day sotorasib covalent capture at Cref
  kcov_adag = 22.0,     # FIT: /day adagrasib covalent capture at Cref
  # --- tumour kinetics ---
  kprol     = 0.026,    # /day max proliferation (at S=1)
  kdeath0   = 0.020,    # /day baseline death  => net untreated growth +0.007/day (doubling ~100 d)
  doce_res  = 0.35,     # fraction of docetaxel kill the resistant clone still receives
  kkill     = 0.032,    # FIT: /day max ERK-suppression-driven death (in addicted cells)
  dep       = 0.33,     # FIT: KRAS-pathway dependence of sensitive cells (oncogene addiction)
  # --- resistance ---
  bypass    = 0.85,     # FIT: de-dependence of resistant clones (near-complete escape)
  fr0       = 0.06,     # FIT: pre-existing resistant fraction of tumour
  kres      = 0.004,    # FIT: /day acquired resistance (sensitive->resistant conversion)
  # --- docetaxel comparator ---
  kkill_doce= 0.009,    # FIT: docetaxel mean daily kill coefficient (per unit normalised exposure)
  kprog_doce= 0.016,    # FIT: /day docetaxel chemoresistance (effect-waning) rate
  kel_doce  = 1.51,     # /day effective elimination (legacy; unused in averaged mode)
  tau_doce  = 21        # day dosing cycle (legacy; unused in averaged mode)
)

# parameters estimated in calibration (the rest are fixed by structure/priors)
FIT_PARS <- c("fb", "kprog", "kcov_soto", "kcov_adag", "kkill", "dep", "bypass", "fr0", "kres", "kkill_doce", "kprog_doce")

# bounds for the fitted parameters
FIT_LB <- c(fb=0.2, kprog=0.002, kcov_soto=2,  kcov_adag=2,  kkill=0.015, dep=0.15, bypass=0.50, fr0=0.02, kres=0.002, kkill_doce=0.002, kprog_doce=0.004)
FIT_UB <- c(fb=4.0, kprog=0.080, kcov_soto=40, kcov_adag=60, kkill=0.12,  dep=0.95, bypass=0.99, fr0=0.30, kres=0.020, kkill_doce=0.060, kprog_doce=0.050)

# ---------------------------------------------------------------------------
# 3. Population variability (for virtual patients). Lognormal multiplicative
#    CVs on selected parameters + a baseline-SLD and PK-exposure factor.
# ---------------------------------------------------------------------------
POP_CV <- c(
  dep    = 0.90,   # KRAS-dependence heterogeneity (logit-scale SD) -- main ORR driver
  kkill  = 0.30,   # response depth heterogeneity
  bypass = 0.25,   # resistance heterogeneity
  fr0    = 0.80,   # pre-existing resistant fraction (highly variable)
  kprol  = 0.30,   # growth-rate heterogeneity
  kcov   = 0.25,   # exposure/engagement heterogeneity (applied to the active drug)
  fb     = 0.30    # feedback strength heterogeneity
)

# ---------------------------------------------------------------------------
# 4. Calibration targets — published summary endpoints (concise, cited).
#    ORR in %, PFS in months. 'w' = weight; 'sd' = approx uncertainty.
#    Sources (see verification/SOURCE_VERIFICATION.md):
#      soto 2L+ NSCLC: CodeBreaK 100 (Skoulidis 2021) ORR 37.1%, mPFS 6.8 mo;
#                      CodeBreaK 200 (de Langen 2023) mPFS 5.6 mo vs docetaxel.
#      adag 2L+ NSCLC: KRYSTAL-1 (Janne 2022) ORR ~43%, mPFS 6.5 mo.
#      docetaxel     : CodeBreaK 200 control ORR ~13%, mPFS 4.5 mo.
#      target engagement: high covalent occupancy at clinical dose (translational).
# ---------------------------------------------------------------------------
TARGETS <- data.table(
  tkey   = c("soto_ORR","soto_PFS","adag_ORR","adag_PFS","doce_ORR","doce_PFS","soto_occ"),
  arm    = c("sotorasib","sotorasib","adagrasib","adagrasib","docetaxel","docetaxel","sotorasib"),
  metric = c("ORR","PFS","ORR","PFS","ORR","PFS","occupancy"),
  value  = c(37.1,      6.8,       43.0,      6.5,       13.2,      4.5,        0.90),
  sd     = c(5.0,       1.0,       6.0,       1.2,       6.0,       1.0,        0.07),
  w      = c(1.0,       1.0,       1.0,       1.0,       0.6,       0.7,        0.6),
  source = c("CodeBreaK100","CodeBreaK100","KRYSTAL-1","KRYSTAL-1","CodeBreaK200","CodeBreaK200","translational")
)

cat("params.R loaded.\n")
