# Data-generating scenarios for the learned-DEF investigation.
# Most scenarios use 2 covariates and are fitted as glm(y ~ x1 + x2); the
# "multivar_*" scenarios use 4 covariates fitted as glm(y ~ x1+x2+x3+x4).
# Each returns list(data, fit).

# alternative families (null handled separately)
ALT_2COV <- c("quad", "quad2", "cubic", "cloglog", "probit",
              "bump", "sine", "step", "interaction", "mixed")
ALT_MULTI <- c("multivar_quad", "multivar_int")
SCENARIOS <- c("null", ALT_2COV, ALT_MULTI)

gen_data <- function(scn, n = 500) {

  # ---- 4-covariate scenarios ----
  if (scn %in% ALT_MULTI) {
    x1 <- runif(n, -3, 3); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- runif(n, -2, 2)
    eta0 <- 0.2 + 0.6 * x1 - 0.4 * x2 + 0.3 * x3 - 0.2 * x4
    p <- switch(scn,
      multivar_quad = plogis(eta0 + 0.22 * x1^2),       # omitted quadratic (4-cov)
      multivar_int  = plogis(eta0 + 0.40 * x1 * x3))     # omitted interaction (4-cov)
    y   <- rbinom(n, 1, p)
    dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
    fit <- suppressWarnings(glm(y ~ x1 + x2 + x3 + x4, data = dat, family = binomial()))
    return(list(data = dat, fit = fit))
  }

  # ---- 2-covariate scenarios ----
  x1 <- runif(n, -3, 3); x2 <- rnorm(n)
  eta0 <- 0.3 + 0.8 * x1 - 0.5 * x2
  p <- switch(scn,
    null        = plogis(eta0),                          # correct model (size)
    quad        = plogis(eta0 + 0.25 * x1^2),            # omitted quadratic in x1
    quad2       = plogis(eta0 + 0.50 * x2^2),            # omitted quadratic in x2
    cubic       = plogis(eta0 + 0.07 * x1^3),            # omitted cubic
    cloglog     = 1 - exp(-exp(0.5 * x1 - 0.3 * x2)),    # wrong link (cloglog)
    probit      = pnorm(0.55 * x1 - 0.35 * x2),          # wrong link (probit)
    bump        = plogis(eta0 + 2.0 * exp(-((x1 - 1)^2) / 0.15)),  # localized bump
    sine        = plogis(eta0 + 1.2 * sin(2 * x1)),      # oscillation
    step        = plogis(eta0 + 1.2 * (x1 > 0.5)),       # threshold / discontinuity
    interaction = plogis(eta0 + 0.4 * x1 * x2),          # omitted interaction
    mixed       = plogis(eta0 + 0.2 * x1^2 + 0.3 * x1 * x2),  # quad + interaction
    ## ---- INDEX-AXIS families, added 2026-08-12 (spec 12) -------------------------
    ## Both are functions of the linear index ALONE: p = f(eta0) with f != plogis.
    ## They exist because the study previously contained NO detectable departure of
    ## this kind -- every index-type family was a link-family swap, and three of those
    ## are undetectable. Amplitudes frozen from the T0 band (spec 12 gate G3) BEFORE
    ## any localization was computed on them: calstrat 1.30 -> global .670,
    ## caltails 0.70 -> see cal_smoke_n1000.csv. Mid-band by the pre-declared rule.
    calstrat    = plogis(eta0 + 1.30 * (plogis(eta0) > 0.75)),   # upper-risk shift
    ## caltails amplitude CORRECTED 1.30->1.00 on 2026-08-12: the value first written (0.70)
    ## was chosen before the T0 sweep landed, and the sweep then showed 0.70 gives global
    ## rejection .175 -- BELOW the .35 band floor declared in spec 12. Of the three amplitudes
    ## swept, only 1.00 is in band (.350). The pre-declared rule therefore selects 1.00; this
    ## is the rule being applied, not a post-hoc choice, and no localization had been computed
    ## on this family when the correction was made.
    caltails    = plogis(eta0 + 1.00 * (plogis(eta0) > 0.80)
                              - 1.00 * (plogis(eta0) < 0.20)),  # both-tail shift
    stop("unknown scenario: ", scn))
  y   <- rbinom(n, 1, p)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- suppressWarnings(glm(y ~ x1 + x2, data = dat, family = binomial()))
  list(data = dat, fit = fit)
}
