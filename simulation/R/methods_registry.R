# =============================================================================
# methods_registry.R -- the method registry.
#
# CHEAP methods (one eval_cheap() call per simulated dataset, all deterministic
# given the data): the fast run.all.gof battery (include_slow = FALSE; all 9
# META_FEATS present), the CORRECTED le-Cessie (v2.4.1, selected explicitly
# because its registry entry is slow=TRUE), and the two covariate-space
# directed bases via the package cdef.gof() (== research methodA_def, verified
# 9/9 identical, Decision D2).
#
# META methods (Stage 2, built from stored per-rep p-values): C.meta (logistic
# on the 9 battery z's), C+A.ridge (glmnet ridge on 11 z's), each with
# held-out-family CV and/or frozen training, calibrated by the conformal
# out-of-sample null reference (and, for the S3 exhibit, the WRONG in-sample
# way as a comparison column).
#
# BAGofT (Zhang, Ding & Yang 2023, JASA) -- headline rival, cost-managed, own
# runner + flag. Published defaults only: nsplits=100, nsim=100, parRF(ntree=60),
# ne = floor(5*sqrt(n)).
# =============================================================================

# ---- cheap-method name sets --------------------------------------------------
META_FEATS9 <- c("HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF",
                 "DEF.poly2", "DEF.poly3", "DEF.stukel")
FEATS11     <- c(META_FEATS9, "A.poly", "A.spline")
CHEAP_METHODS <- c(META_FEATS9, "Stukel", "Osius-Rojek", "le-Cessie",
                   "CCT", "A.poly", "A.spline")
# CCT = run.all.gof's "Ensemble.Vote(3DEF)" (Cauchy combination, Liu-Xie 2020).

# ---- one dataset -> all cheap p-values --------------------------------------
eval_cheap <- function(fit) {
  out <- setNames(rep(NA_real_, length(CHEAP_METHODS)), CHEAP_METHODS)
  r <- tryCatch(suppressWarnings(run.all.gof(fit, include_slow = FALSE, install = "no")),
                error = function(e) NULL)
  if (!is.null(r)) {
    pv <- setNames(r$p_value, r$Test)
    grab <- c(META_FEATS9, "Stukel", "Osius-Rojek")
    out[grab] <- pv[grab]
    out["CCT"] <- unname(pv["Ensemble.Vote(3DEF)"])
  }
  # corrected le-Cessie (v2.4.1 tree); slow-flagged in the registry, so select it
  out["le-Cessie"] <- tryCatch(
    suppressWarnings(suppressMessages(
      run.all.gof(fit, tests = "le-Cessie", include_slow = TRUE, install = "no")$p_value)),
    error = function(e) NA_real_)
  out["A.poly"]   <- tryCatch(suppressWarnings(cdef.gof(fit, basis = "poly")$p_value),
                              error = function(e) NA_real_)
  out["A.spline"] <- tryCatch(suppressWarnings(cdef.gof(fit, basis = "spline")$p_value),
                              error = function(e) NA_real_)
  out
}

# ---- evidence transform + meta-learners (ported from run_confirmatory.R) ----
z_evidence <- function(p) {           # one-sided evidence, larger = more misfit
  z <- qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6), lower.tail = FALSE)
  z[!is.finite(z) | is.na(z)] <- 0
  z
}
feats_matrix <- function(pmat, which_feats) {
  # pmat: data.frame/matrix of per-rep p-values with CHEAP_METHODS columns
  Z <- apply(as.matrix(pmat[, which_feats, drop = FALSE]), 2, z_evidence)
  colnames(Z) <- make.names(which_feats)
  Z
}
emp_p <- function(scores, ref)        # conformal one-sided empirical p
  (1 + colSums(outer(ref, scores, ">="))) / (length(ref) + 1)

train_logit <- function(FM, y) {
  d <- as.data.frame(FM); d$.lab <- y
  suppressWarnings(glm(.lab ~ ., data = d, family = binomial()))
}
score_logit <- function(m, FM)
  as.numeric(predict(m, newdata = as.data.frame(FM), type = "link"))
train_ridge <- function(FM, y, nfolds = 10L) {
  nfolds <- max(3L, min(nfolds, floor(length(y) / 3)))
  glmnet::cv.glmnet(as.matrix(FM), y, family = "binomial", alpha = 0, nfolds = nfolds)
}
score_ridge <- function(m, FM)
  as.numeric(predict(m, newx = as.matrix(FM), s = "lambda.min", type = "link"))

# ---- BAGofT wrapper (published defaults; guarded by JASA$RUN_BAGOFT) --------
# dat is generated OUTSIDE (identical per-cell seeds as the cheap grid);
# BAGofT's internal RNG is re-seeded by the caller with METHOD_RNG_OFFSET.
# NOTE (KB sec. 7): pass a FRESH data.frame, never model.frame(fit) -- a stale
# `terms` attribute makes BAGofT return p = 0.
bagoft_p <- function(dat, formula = y ~ x1 + x2,
                     nsplits = JASA$BAG_NSPLITS, nsim = JASA$BAG_NSIM,
                     ntree = JASA$BAG_NTREE) {
  environment(formula) <- globalenv()
  dat <- as.data.frame(lapply(dat, as.vector))   # strip any attributes
  t0 <- Sys.time()
  b <- tryCatch({
    invisible(capture.output(
      res <- suppressMessages(suppressWarnings(
        BAGofT::BAGofT(testModel = BAGofT::testGlmBi(formula = formula, link = "logit"),
                       parFun = BAGofT::parRF(ntree = ntree),
                       data = dat,
                       nsplits = nsplits,
                       ne = floor(5 * sqrt(nrow(dat))),
                       nsim = nsim)))))
    res
  }, error = function(e) NULL)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(b)) return(list(p = NA_real_, p_med = NA_real_, p_min = NA_real_,
                              stat_min = NA_real_, time_s = el))
  # BAGofT returns three bootstrap p-values -- p.value (splits combined by the mean, the
  # paper's headline), p.value2 (median), p.value3 (minimum) -- and, separately, the combined
  # STATISTICS pmean / pmedian / pmin. Before 2026-09-16 this function returned pmin as
  # p_min. pmin is the smallest per-split p-value: it lies in [0,1] and looks like a p-value,
  # but over 100 splits it is almost always below .05. It is kept, correctly named, as stat_min.
  getnum <- function(x) if (!is.null(x) && is.numeric(x) && length(x) == 1) as.numeric(x) else NA_real_
  list(p        = getnum(b$p.value),
       p_med    = getnum(b$p.value2),
       p_min    = getnum(b$p.value3),
       stat_min = getnum(b$pmin),
       time_s   = el)
}
