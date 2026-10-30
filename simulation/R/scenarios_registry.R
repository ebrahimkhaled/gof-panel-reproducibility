# =============================================================================
# scenarios_registry.R -- the pre-registered family sets (Evidence Plan sec. 2)
# and one unified generator gen_scenario(scn, n) -> list(data, fit, formula).
#
# The 12 original families come UNCHANGED from learned-def/R/scenarios.R
# (sourced, not copied -- guarantees identical DGPs). The 5 richer families are
# ported verbatim from exp_richer_families.R. The 3 deployment families (D1-D3)
# and the t8link all-fail family are NEW here, exactly as pre-registered in
# 02_EVIDENCE_PLAN.md section 2.
# =============================================================================

# ---- family sets (pre-registered) -------------------------------------------
POOL_OLD12  <- c("quad", "quad2", "cubic", "cloglog", "probit", "bump", "sine",
                 "step", "interaction", "mixed", "multivar_quad", "multivar_int")
POOL_NEW5   <- c("logabs", "threshx2", "dbump", "skewlink", "quadint_outlier")
## INDEX-axis families (spec 12, 2026-08-12): p = f(eta0), detectable, NOT link-family swaps.
## Evaluation-only -- deliberately NOT in POOL_P, so the frozen rule never trains on them.
POOL_INDEX2 <- c("calstrat", "caltails")
POOL_P      <- c(POOL_OLD12, POOL_NEW5)          # training pool (17 alts)
DEPLOY_D    <- c("cauchit", "omitvar", "etatrans")  # frozen-combiner deployment set
ALLFAIL_F   <- c("probit", "cubic", "t8link")       # mechanism block (+ threshx2 borderline)
ALL_ALT     <- c(POOL_P, DEPLOY_D, "t8link", POOL_INDEX2)
ALL_FAMILIES <- c("null", ALL_ALT)
CONFIRM_ALT <- POOL_OLD12                        # S2 confirmatory 12-family set

# t8link latent scale: s such that Var(s * t_8) = pi^2/3 (the logistic latent
# variance): s = sqrt((pi^2/3) / (8/6)) = pi/2. Deliberately near-logistic =>
# near-invisible (all-fail by design). FLAGGED for T0.5 amplitude confirmation.
T8_SCALE <- pi / 2

# amplitude knobs for the deployment set (T0.5 smoke-tuning adjusts these ONLY)
D_AMP <- list(cauchit  = c(a1 = 0.45, a2 = -0.30),
              omitvar  = c(sd_u = 0.8),
              etatrans = c(pow = 1.3),
              t8link   = c(a1 = 0.55, a2 = -0.35))

# ---- the 5 richer families (verbatim port of exp_richer_families.R) ---------
gen_new <- function(scn, n = 500) {
  x1 <- runif(n, -3, 3); x2 <- rnorm(n)
  eta0 <- 0.3 + 0.8 * x1 - 0.5 * x2
  p <- switch(scn,
    logabs          = plogis(eta0 + 0.9 * log(abs(x1) + 0.5)),
    threshx2        = plogis(eta0 + 1.1 * (x2 > 0)),
    dbump           = plogis(eta0 + 1.8 * exp(-((x1 + 1.3)^2) / 0.20)
                                  - 1.8 * exp(-((x1 - 1.3)^2) / 0.20)),
    skewlink        = 1 - exp(-exp(0.45 * x1 - 0.30 * x2 + 0.18 * x1^2 * (x1 < 0))),
    quadint_outlier = {
      k <- rbinom(n, 1, 0.05)
      x1 <- ifelse(k == 1, rt(n, df = 2) * 2, x1)
      eta0 <- 0.3 + 0.8 * x1 - 0.5 * x2
      plogis(eta0 + 0.22 * x1^2 + 0.35 * x1 * x2)
    },
    stop("unknown richer scenario: ", scn))
  y   <- rbinom(n, 1, p)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- suppressWarnings(glm(y ~ x1 + x2, data = dat, family = binomial()))
  list(data = dat, fit = fit)
}

# ---- deployment set D1-D3 + t8link (NEW; pre-registered definitions) --------
gen_deploy <- function(scn, n = 500) {
  x1 <- runif(n, -3, 3); x2 <- rnorm(n)
  eta0 <- 0.3 + 0.8 * x1 - 0.5 * x2
  p <- switch(scn,
    # D1: cauchit link (heavy-tailed link misspecification)
    cauchit  = pcauchy(D_AMP$cauchit["a1"] * x1 + D_AMP$cauchit["a2"] * x2),
    # D2: omitted independent covariate (unobserved heterogeneity; hard by design)
    omitvar  = plogis(eta0 + D_AMP$omitvar["sd_u"] * rnorm(n)),
    # D3: monotone eta transformation
    etatrans = plogis(sign(eta0) * abs(eta0)^D_AMP$etatrans["pow"]),
    # F: latent-t(8) link, scaled near-logistic (all-fail by design)
    t8link   = pt((D_AMP$t8link["a1"] * x1 + D_AMP$t8link["a2"] * x2) / T8_SCALE, df = 8),
    stop("unknown deployment scenario: ", scn))
  p   <- pmin(pmax(as.numeric(p), 1e-12), 1 - 1e-12)
  y   <- rbinom(n, 1, p)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- suppressWarnings(glm(y ~ x1 + x2, data = dat, family = binomial()))
  list(data = dat, fit = fit)
}

# ---- unified generator -------------------------------------------------------
# Requires learned-def/R/scenarios.R (gen_data, ALT_2COV, ALT_MULTI) sourced.
gen_scenario <- function(scn, n = 500) {
  g <- if (scn %in% POOL_INDEX2)                  gen_data(scn, n)
       else if (scn %in% POOL_NEW5)                    gen_new(scn, n)
       else if (scn %in% c(DEPLOY_D, "t8link"))   gen_deploy(scn, n)
       else                                       gen_data(scn, n)   # null + old 12
  g$formula <- if (scn %in% ALT_MULTI) y ~ x1 + x2 + x3 + x4 else y ~ x1 + x2
  g
}
