## =============================================================================
## select_final_heldout.R -- FINAL SELECT PHASE (2026-08-07)
##
## ONE rule, written down in results/jasa_stage2/SELECTION_PREREG.md BEFORE
## this script was written, is evaluated ONCE on the held-out test_* families.
##
## SELECTED RULE: NNCauchy (nested-CV, select=mean)
##   score  S = sum_j w_j * tan((0.5 - clip(p_j)) * pi),  j in FEATS11
##   w_j    = ahat_j^gamma / sum_k ahat_k^gamma   (w >= 0, sum w = 1)
##   ahat_j = pmax( agg_f HIT[f,j], 0.005 ),  HIT = per-test size-adj marginal
##            power on training family f (cutoff = 95th pct of conformal null)
##   (gamma, agg) chosen by INNER leave-one-TRAINING-family-out CV, objective =
##            mean size-adjusted power, ties -> worst-case -> smaller gamma
##   gamma in {0,.25,.5,1,2,3,5,Inf}; agg in {mean,min}; (0,min) dropped
##   gamma = 0 reproduces CCT-full EXACTLY.
##
## Fold structure = run_stage2_meta.R section 1, verbatim:
##   for f in CONFIRM_ALT: fit on train_null + train_{CONFIRM_ALT \ f},
##   score test_f and test_null, conformal p from conformal_null via emp_p(),
##   size = mean(p_null <= .05), size-adj power = mean(p_f <= q_.05(p_null)).
##
## Comparators recomputed on the SAME held-out families: CCT-full (equal weight)
## and C+A.ridge(HO) (the shipped learned rule).
##
## No new simulation. Deterministic apart from cv.glmnet inside C+A.ridge, which
## is seeded.
##
## Usage:  Rscript select_final_heldout.R
## Output: results/jasa_stage2/FINAL_heldout_n{500,1000}.csv
##         results/jasa_stage2/FINAL_heldout_selection.csv
##         results/jasa_stage2/FINAL_heldout_weights.csv
##         results/jasa_stage2/FINAL_heldout_LOG.txt
## =============================================================================

suppressWarnings(suppressMessages(library(glmnet)))

BASE <- "."
S1   <- file.path(BASE, "results", "jasa_stage1")
OUT  <- file.path(BASE, "results", "jasa_stage2")
LOG  <- file.path(OUT, "FINAL_heldout_LOG.txt")
cat("", file = LOG)
say <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = LOG, append = TRUE) }

ALPHA <- 0.05
NS    <- c(500L, 1000L)
SEED  <- 20260807L

META_FEATS9 <- c("HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF",
                 "DEF.poly2", "DEF.poly3", "DEF.stukel")
FEATS11     <- c(META_FEATS9, "A.poly", "A.spline")
CHEAP_METHODS <- c(META_FEATS9, "Stukel", "Osius-Rojek", "le-Cessie",
                   "CCT", "A.poly", "A.spline")
POOL_OLD12  <- c("quad", "quad2", "cubic", "cloglog", "probit", "bump", "sine",
                 "step", "interaction", "mixed", "multivar_quad", "multivar_int")
CONFIRM_ALT <- POOL_OLD12
ALT11       <- setdiff(POOL_OLD12, "probit")
K           <- length(FEATS11)

## ---- transforms, verbatim from the shipped code -----------------------------
z_evidence <- function(p) {                       # methods_registry.R
  z <- qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6), lower.tail = FALSE)
  z[!is.finite(z) | is.na(z)] <- 0
  z
}
clipp  <- function(p) pmin(pmax(p, 1e-12), 1 - 1e-12)
tan_of <- function(p) {                           # run_stage2_meta.R cct_full_p
  t <- tan((0.5 - clipp(p)) * pi)
  t[!is.finite(t) | is.na(t)] <- 0                # NA -> neutral, mirrors z_evidence
  t
}
emp_p <- function(scores, ref)                    # methods_registry.R
  (1 + colSums(outer(ref, scores, ">="))) / (length(ref) + 1)

read_cell <- function(role, fam, n) {
  f <- file.path(S1, sprintf("%s_%s_n%d_pvalues.csv", role, fam, n))
  if (!file.exists(f)) stop("missing cell: ", f)
  utils::read.csv(f, check.names = FALSE)
}
prep <- function(d) {
  P <- as.matrix(d[, FEATS11, drop = FALSE])
  Z <- apply(P, 2, z_evidence); Tt <- apply(P, 2, tan_of)
  colnames(Z) <- colnames(Tt) <- make.names(FEATS11)
  list(Z = Z, T = Tt, nrep = nrow(P), nNA = sum(!is.finite(P) | is.na(P)))
}
pow_at <- function(s_alt, s_null, alpha = ALPHA)
  mean(s_alt > stats::quantile(s_null, 1 - alpha, na.rm = TRUE, names = FALSE), na.rm = TRUE)

## =============================================================================
all_tabs <- list(); all_sel <- list(); all_w <- list()

for (n in NS) {

  say("\n=============================== n = %d ===============================", n)

  CAL <- prep(read_cell("conformal", "null", n))
  TR  <- lapply(setNames(c("null", POOL_OLD12), c("null", POOL_OLD12)),
                function(f) prep(read_cell("train", f, n)))
  TEraw <- lapply(setNames(c("null", POOL_OLD12), c("null", POOL_OLD12)),
                  function(f) read_cell("test", f, n))
  TE  <- lapply(TEraw, prep)

  say("cells: conformal null %d | train null %d, 12 train alts x %d | test null %d, test alts x %d",
      CAL$nrep, TR$null$nrep, TR$quad$nrep, TE$null$nrep, TE$quad$nrep)
  say("non-finite FEATS11 entries: conformal %d, train(all) %d, test(all) %d",
      CAL$nNA, sum(vapply(TR, `[[`, numeric(1), "nNA")), sum(vapply(TE, `[[`, numeric(1), "nNA")))

  ## ---- ahat machinery (training families + conformal null ONLY) -------------
  cthr <- apply(CAL$T, 2, stats::quantile, probs = 1 - ALPHA, names = FALSE)
  HIT  <- t(vapply(POOL_OLD12, function(f)
             colMeans(sweep(TR[[f]]$T, 2, cthr, ">")), numeric(K)))
  rownames(HIT) <- POOL_OLD12; colnames(HIT) <- make.names(FEATS11)

  FLOOR <- 0.005
  ahat_of <- function(fams, agg)
    pmax(if (agg == "mean") colMeans(HIT[fams, , drop = FALSE])
         else               apply(HIT[fams, , drop = FALSE], 2, min), FLOOR)
  w_of <- function(a, gamma) {
    if (gamma == 0) return(rep(1 / K, K))
    if (!is.finite(gamma)) { w <- as.numeric(a == max(a)); return(w / sum(w)) }
    w <- a^gamma; w / sum(w)
  }
  GAMMAS <- c(0, 0.25, 0.5, 1, 2, 3, 5, Inf)
  CFG <- expand.grid(gamma = GAMMAS, agg = c("mean", "min"), stringsAsFactors = FALSE)
  CFG <- CFG[!(CFG$gamma == 0 & CFG$agg == "min"), ]

  ## ==========================================================================
  ## OUTER: held-out-family CV over CONFIRM_ALT, evaluated on test_*
  ## ==========================================================================
  sz  <- list(); pw <- list()
  for (rl in c("NNCauchy", "CCT-full(conformal)", "C+A.ridge(HO)")) {
    sz[[rl]] <- setNames(numeric(length(CONFIRM_ALT)), CONFIRM_ALT)
    pw[[rl]] <- setNames(numeric(length(CONFIRM_ALT)), CONFIRM_ALT)
  }

  for (f in CONFIRM_ALT) {
    trf <- setdiff(CONFIRM_ALT, f)

    ## ---------- SELECTED RULE: inner LOTFO over the in-fold families --------
    inner <- matrix(NA_real_, nrow = nrow(CFG), ncol = length(trf))
    for (gi in seq_along(trf)) {
      g_out <- trf[gi]; S_in <- setdiff(trf, g_out)
      a_mean <- ahat_of(S_in, "mean"); a_min <- ahat_of(S_in, "min")
      for (ci in seq_len(nrow(CFG))) {
        a <- if (CFG$agg[ci] == "mean") a_mean else a_min
        w <- w_of(a, CFG$gamma[ci])
        inner[ci, gi] <- pow_at(as.numeric(TR[[g_out]]$T %*% w),
                                as.numeric(CAL$T %*% w))
      }
    }
    im <- rowMeans(inner); iw <- apply(inner, 1, min)
    best <- order(-round(im, 6), -round(iw, 6), CFG$gamma)[1]
    a <- ahat_of(trf, CFG$agg[best])
    w <- w_of(a, CFG$gamma[best])

    nref <- as.numeric(CAL$T %*% w)
    p0   <- emp_p(as.numeric(TE$null$T %*% w), nref)
    pf   <- emp_p(as.numeric(TE[[f]]$T %*% w), nref)
    sz[["NNCauchy"]][f] <- mean(p0 <= ALPHA, na.rm = TRUE)
    pw[["NNCauchy"]][f] <- mean(pf <= quantile(p0, ALPHA, na.rm = TRUE), na.rm = TRUE)

    all_sel[[length(all_sel) + 1L]] <- data.frame(
      n = n, fold = f, gamma = CFG$gamma[best], agg = CFG$agg[best],
      inner_mean = im[best], inner_worst = iw[best],
      w_min = min(w), w_max = max(w), w_max_feature = make.names(FEATS11)[which.max(w)])
    all_w[[length(all_w) + 1L]] <- data.frame(
      n = n, fold = f, feature = make.names(FEATS11), w = as.numeric(w),
      ahat = as.numeric(a))

    ## ---------- COMPARATOR 1: CCT-full, equal weight, conformal ------------
    we <- rep(1 / K, K)
    nrefe <- as.numeric(CAL$T %*% we)
    p0e <- emp_p(as.numeric(TE$null$T %*% we), nrefe)
    pfe <- emp_p(as.numeric(TE[[f]]$T %*% we), nrefe)
    sz[["CCT-full(conformal)"]][f] <- mean(p0e <= ALPHA, na.rm = TRUE)
    pw[["CCT-full(conformal)"]][f] <- mean(pfe <= quantile(p0e, ALPHA, na.rm = TRUE), na.rm = TRUE)

    ## ---------- COMPARATOR 2: C+A.ridge(HO), the shipped learned rule ------
    Xtr <- do.call(rbind, c(list(TR$null$Z), lapply(trf, function(g) TR[[g]]$Z)))
    ytr <- c(rep(0L, nrow(TR$null$Z)),
             unlist(lapply(trf, function(g) rep(1L, nrow(TR[[g]]$Z)))))
    set.seed(SEED)
    nf <- max(3L, min(10L, floor(length(ytr) / 3)))
    m  <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 0, nfolds = nf)
    scr <- function(M) as.numeric(predict(m, newx = M, s = "lambda.min", type = "link"))
    nrefr <- scr(CAL$Z)
    p0r <- emp_p(scr(TE$null$Z), nrefr); pfr <- emp_p(scr(TE[[f]]$Z), nrefr)
    sz[["C+A.ridge(HO)"]][f] <- mean(p0r <= ALPHA, na.rm = TRUE)
    pw[["C+A.ridge(HO)"]][f] <- mean(pfr <= quantile(p0r, ALPHA, na.rm = TRUE), na.rm = TRUE)
  }

  ## ---- CCT-full the OFFICIAL way (classical row: raw p, test-null cutoff) --
  cct_p <- function(d) {
    Tst <- rowMeans(sapply(FEATS11, function(m) tan_of(d[[m]])))
    0.5 - atan(Tst) / pi
  }
  pn_cct <- cct_p(TEraw$null)
  cut_cct <- quantile(pn_cct, ALPHA, na.rm = TRUE, names = FALSE)
  cct_off <- vapply(CONFIRM_ALT, function(a) mean(cct_p(TEraw[[a]]) <= cut_cct, na.rm = TRUE), numeric(1))
  size_cct_off <- mean(pn_cct < ALPHA, na.rm = TRUE)

  ## ---- ORACLE: per-family hindsight best single deposited method -----------
  pn <- TEraw$null
  cut_c <- vapply(CHEAP_METHODS, function(m)
    quantile(pn[[m]], ALPHA, na.rm = TRUE, names = FALSE), numeric(1))
  powc <- t(vapply(CHEAP_METHODS, function(m)
    vapply(CONFIRM_ALT, function(a) mean(TEraw[[a]][[m]] <= cut_c[m], na.rm = TRUE),
           numeric(1)), numeric(length(CONFIRM_ALT))))
  oracle <- apply(powc, 2, max)
  oracle_mean <- mean(oracle[ALT11])

  ## ---- assemble ------------------------------------------------------------
  mk <- function(name, size, p) data.frame(
    Method = name, n = n, Size = round(size, 4),
    as.list(round(p[CONFIRM_ALT], 4)),
    mean_ALT11 = round(mean(p[ALT11]), 4), worst_ALT11 = round(min(p[ALT11]), 4),
    oracle_gap = round(oracle_mean - mean(p[ALT11]), 4), check.names = FALSE)

  tab <- rbind(
    mk("NNCauchy (SELECTED)",            mean(sz[["NNCauchy"]]),            pw[["NNCauchy"]]),
    mk("CCT-full (equal wt, conformal)", mean(sz[["CCT-full(conformal)"]]), pw[["CCT-full(conformal)"]]),
    mk("CCT-full (official classical row)", size_cct_off,                   cct_off),
    mk("C+A.ridge(HO) (shipped learned)", mean(sz[["C+A.ridge(HO)"]]),      pw[["C+A.ridge(HO)"]]),
    mk("ORACLE (hindsight best single)", NA_real_,                          oracle))
  rownames(tab) <- NULL
  utils::write.csv(tab, file.path(OUT, sprintf("FINAL_heldout_n%d.csv", n)), row.names = FALSE)
  all_tabs[[as.character(n)]] <- tab

  say("\n--- HELD-OUT (test_*) RESULTS, n = %d, alpha = %.2f, probit excluded ---", n, ALPHA)
  print(tab[, c("Method", "Size", "mean_ALT11", "worst_ALT11", "oracle_gap")])
  capture.output(print(tab[, c("Method", "Size", "mean_ALT11", "worst_ALT11", "oracle_gap")]),
                 file = LOG, append = TRUE)
  say("\nper-family:")
  capture.output(print(tab[, c("Method", CONFIRM_ALT)]), file = LOG, append = TRUE)
  print(tab[, c("Method", CONFIRM_ALT)])

  ## ---- paired inference vs CCT-full over the 11 held-out families ----------
  d1 <- pw[["NNCauchy"]][ALT11] - pw[["CCT-full(conformal)"]][ALT11]
  d2 <- pw[["NNCauchy"]][ALT11] - pw[["C+A.ridge(HO)"]][ALT11]
  fam_t <- function(d) { m <- mean(d); s <- sd(d) / sqrt(length(d))
    c(mean = m, se = s, t = if (s > 0) m / s else NA_real_,
      up = sum(d > 0), down = sum(d < 0), tie = sum(d == 0)) }
  say("\npaired over the 11 held-out families:")
  say("  NNCauchy - CCT-full : mean %+0.4f  se %0.4f  t %+0.2f  up %d / down %d / tie %d",
      fam_t(d1)[1], fam_t(d1)[2], fam_t(d1)[3], fam_t(d1)[4], fam_t(d1)[5], fam_t(d1)[6])
  say("  NNCauchy - C+A.ridge: mean %+0.4f  se %0.4f  t %+0.2f  up %d / down %d / tie %d",
      fam_t(d2)[1], fam_t(d2)[2], fam_t(d2)[3], fam_t(d2)[4], fam_t(d2)[5], fam_t(d2)[6])
  nrep_f <- TE$quad$nrep
  say("  per-family MC se of a single power estimate ~ %0.4f (%d reps)",
      sqrt(0.25 / nrep_f), nrep_f)

  ## ---- certificate: gamma = 0 reproduces CCT-full exactly ------------------
  w0 <- rep(1 / K, K)
  chk <- max(abs(pw[["CCT-full(conformal)"]][CONFIRM_ALT] -
    vapply(CONFIRM_ALT, function(f) {
      nr <- as.numeric(CAL$T %*% w0)
      p0 <- emp_p(as.numeric(TE$null$T %*% w0), nr)
      mean(emp_p(as.numeric(TE[[f]]$T %*% w0), nr) <= quantile(p0, ALPHA), na.rm = TRUE)
    }, numeric(1))))
  say("  CERTIFICATE gamma=0 == CCT-full : max |power diff| = %.3e", chk)
  say("  official-row vs conformal CCT-full: max |power diff| = %.4f",
      max(abs(cct_off[CONFIRM_ALT] - pw[["CCT-full(conformal)"]][CONFIRM_ALT])))
}

utils::write.csv(do.call(rbind, all_sel), file.path(OUT, "FINAL_heldout_selection.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, all_w),   file.path(OUT, "FINAL_heldout_weights.csv"),   row.names = FALSE)

say("\n--- selected (gamma, agg) per production fold ---")
capture.output(print(do.call(rbind, all_sel)[, c("n","fold","gamma","agg","w_max_feature","w_max")]),
               file = LOG, append = TRUE)
print(do.call(rbind, all_sel)[, c("n","fold","gamma","agg","w_max_feature","w_max")])
say("\nDONE.")
