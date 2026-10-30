## =============================================================================
## cand_nonneg.R -- CANDIDATE: NON-NEGATIVE-WEIGHT STACKING  (2026-08-07)
##
## MOTIVATION.  Every coordinate the learner sees is oriented the same way
## (higher = more misfit).  A combiner that puts a NEGATIVE weight on one of
## them is asserting "this test firing is evidence of FIT", which no member of
## the battery supports.  Diagnostic A showed the shipped ridge does exactly
## that on 2-3 coordinates per fold (HL, Pigeon-Heyse at n=500; DEF.poly2 at
## n=1000), with |sd/mean| up to 0.87 -- i.e. sign-unstable noise-fitting that
## can cancel genuine signal.  This script constrains beta >= 0.
##
## TWO COORDINATE SYSTEMS ARE TRIED, because non-negativity alone is not the
## whole story:
##
##   (P) PROBIT coordinates  z_j = Phi^{-1}(1 - p_j), clip 1e-6 (the SHIPPED map).
##       Rule: score = sum_j w_j z_j, w >= 0, fitted by glmnet(alpha = 0,
##       lower.limits = 0) -- non-negative ridge logistic stacking (NNLS-style).
##       This is the literal brief.  Its equal-weight member is Zmean, NOT
##       CCT-full, so the class does NOT contain the baseline.
##
##   (C) CAUCHY coordinates  t_j = tan((0.5 - p_j) * pi), clip 1e-12 (exactly
##       the map the CCT-full baseline aggregates).
##       Rule: score = sum_j w_j t_j with w >= 0 and sum_j w_j = 1, i.e. a
##       PROPER WEIGHTED AVERAGE directly comparable to equal weights.
##       *** w = (1/K, ..., 1/K) reproduces CCT-full EXACTLY. ***
##       So this non-negative class CONTAINS the baseline by construction and a
##       correctly-selected member cannot be beaten by it except through
##       selection noise.  This is the only version of "non-negative stacking"
##       with a certificate against CCT-full.
##
##       The weights are the one-parameter detectability tilt
##             w_j  propto  ahat_j ^ gamma ,        gamma >= 0,
##       where ahat_j is the training-set marginal size-adjusted power of test j
##       (aggregated over training families by mean or by min).  gamma = 0 IS
##       CCT-full; gamma -> large is "use only the most detectable test".
##       gamma and the aggregator are chosen by an INNER leave-one-training-
##       family-out CV nested inside every outer fold, so the reported number
##       estimates the whole PROCEDURE, not a hand-picked hyper-parameter.
##
## HONEST PROTOCOL (binding).  No test_* cell is opened anywhere in this file.
##   * weight fitting / hyper-parameter selection -> train_* families only
##   * evaluation -> HELD-OUT TRAINING family (LOTFO), never test_*
##   * null calibration -> conformal_null_* (null-only; carries no information
##     about which alternative is which, and is the declared calibration set
##     available at training time)
## The outer fold structure mirrors the shipped rule in run_stage2_meta.R
## (train on null + POOL_OLD12 \ {f}); only the EVALUATION cell is swapped from
## test_f to train_f.  The CCT-full baseline is recomputed HERE on the SAME
## folds, the SAME complete-case rows and the SAME calibration set, so the
## comparison is like-for-like.
##
## Scoring convention: higher = more misfit; size-adjusted power =
## P(score > (1 - alpha) quantile of the conformal-null scores), alpha = 0.05.
##
## No new simulation: pure post-processing of the deposited per-rep p-values.
##
## Usage:  Rscript cand_nonneg.R
## Output: results/jasa_stage2/cand_nonneg_trainCV_n<n>.csv
##         results/jasa_stage2/cand_nonneg_weights_n<n>.csv
##         results/jasa_stage2/cand_nonneg_log.txt
## =============================================================================

suppressWarnings(suppressMessages(library(glmnet)))

BASE <- "."
S1   <- file.path(BASE, "results", "jasa_stage1")
OUT  <- file.path(BASE, "results", "jasa_stage2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG  <- file.path(OUT, "cand_nonneg_log.txt")
cat("", file = LOG)
say <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = LOG, append = TRUE) }

ALPHA <- 0.05
NS    <- c(500L, 1000L)
SEED  <- 20260807L              # same seed as diagA, so ridge rows reproduce

META_FEATS9 <- c("HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF",
                 "DEF.poly2", "DEF.poly3", "DEF.stukel")
FEATS11     <- c(META_FEATS9, "A.poly", "A.spline")
POOL_OLD12  <- c("quad", "quad2", "cubic", "cloglog", "probit", "bump", "sine",
                 "step", "interaction", "mixed", "multivar_quad", "multivar_int")
ALT11       <- setdiff(POOL_OLD12, "probit")   # probit excluded from summaries
K           <- length(FEATS11)

## ---- exact reproductions of the shipped transforms --------------------------
z_evidence <- function(p) {                      # methods_registry.R, clip 1e-6
  z <- qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6), lower.tail = FALSE)
  z[!is.finite(z) | is.na(z)] <- 0
  z
}
clipp  <- function(p) pmin(pmax(p, 1e-12), 1 - 1e-12)
tan_of <- function(p) tan((0.5 - clipp(p)) * pi)  # run_stage2_meta.R cct_full_p

read_cell <- function(role, fam, n) {
  f <- file.path(S1, sprintf("%s_%s_n%d_pvalues.csv", role, fam, n))
  if (!file.exists(f)) return(NULL)
  utils::read.csv(f, check.names = FALSE)
}
## one cell -> the two coordinate matrices on the SAME complete-case rows
prep <- function(d) {
  P <- as.matrix(d[, FEATS11, drop = FALSE])
  keep <- stats::complete.cases(P) & apply(is.finite(P), 1, all)
  P <- P[keep, , drop = FALSE]
  Z <- apply(P, 2, z_evidence); T <- apply(P, 2, tan_of)
  colnames(Z) <- colnames(T) <- make.names(FEATS11)
  list(Z = Z, T = T, nrep = nrow(P))
}
pow_at <- function(s_alt, s_null, alpha = ALPHA)
  mean(s_alt > stats::quantile(s_null, 1 - alpha, na.rm = TRUE, names = FALSE), na.rm = TRUE)

## =============================================================================
for (n in NS) {

  say("\n================================ n = %d ================================", n)
  CALd <- read_cell("conformal", "null", n); if (is.null(CALd)) { say("no conformal null"); next }
  CAL  <- prep(CALd)
  TRd  <- lapply(setNames(c("null", POOL_OLD12), c("null", POOL_OLD12)),
                 read_cell, role = "train", n = n)
  stopifnot(all(vapply(TRd, Negate(is.null), logical(1))))
  TR <- lapply(TRd, prep)
  say("cells: train null = %d reps, %d training alt families x %d reps, conformal null = %d reps",
      TR$null$nrep, length(POOL_OLD12), TR$quad$nrep, CAL$nrep)

  ## ==========================================================================
  ## PRE-COMPUTE: per-test size-adjusted marginal power on every TRAINING family
  ## (thresholds from the conformal null; this is ahat's raw material)
  ## ==========================================================================
  cthr <- apply(CAL$T, 2, stats::quantile, probs = 1 - ALPHA, names = FALSE)  # per-test 95% null quantile
  HIT  <- t(vapply(POOL_OLD12, function(f)
             colMeans(sweep(TR[[f]]$T, 2, cthr, ">")), numeric(K)))
  rownames(HIT) <- POOL_OLD12; colnames(HIT) <- make.names(FEATS11)

  ## ---- the tilt map: ahat -> simplex weights -------------------------------
  FLOOR <- 0.005                                  # keeps ahat^gamma well defined
  ahat_of <- function(fams, agg)
    pmax(if (agg == "mean") colMeans(HIT[fams, , drop = FALSE])
         else               apply(HIT[fams, , drop = FALSE], 2, min), FLOOR)
  w_of <- function(a, gamma) {
    if (gamma == 0) return(rep(1 / K, K))         # <-- EXACTLY CCT-full
    if (!is.finite(gamma)) { w <- as.numeric(a == max(a)); return(w / sum(w)) }
    w <- a^gamma; w / sum(w)
  }
  GAMMAS <- c(0, 0.25, 0.5, 1, 2, 3, 5, Inf)
  AGGS   <- c("mean", "min")
  CFG    <- expand.grid(gamma = GAMMAS, agg = AGGS, stringsAsFactors = FALSE)
  CFG    <- CFG[!(CFG$gamma == 0 & CFG$agg == "min"), ]   # gamma=0 identical for both aggs

  ## ---- crude alternative: EQUAL weights over the top-m tests by ahat --------
  ## (probe: is any gain really just SUBSET SELECTION rather than weighting?)
  w_topm <- function(a, m) { w <- as.numeric(rank(-a, ties.method = "first") <= m); w / sum(w) }
  MS <- seq_len(K)

  ## ==========================================================================
  ## OUTER LOOP: leave-one-TRAINING-family-out
  ## ==========================================================================
  rows <- list(); wrows <- list(); sel <- list()
  for (f in POOL_OLD12) {
    trf <- setdiff(POOL_OLD12, f)
    Zf  <- TR[[f]]$Z; Tf <- TR[[f]]$T

    ## ---------------- (P) probit-coordinate stacking, glmnet ---------------
    X  <- do.call(rbind, c(list(TR$null$Z), lapply(trf, function(g) TR[[g]]$Z)))
    y  <- c(rep(0L, nrow(TR$null$Z)),
            unlist(lapply(trf, function(g) rep(1L, nrow(TR[[g]]$Z)))))
    set.seed(SEED)
    m_un <- cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 10)
    set.seed(SEED)
    m_nn <- cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 10,
                      lower.limits = 0)
    ## essentially-unpenalised non-negative logistic (end of the glmnet path)
    g_nn <- glmnet(X, y, family = "binomial", alpha = 0, lower.limits = 0,
                   lambda = c(1, 1e-2, 1e-4, 1e-6))
    lam_end <- min(g_nn$lambda)

    sc <- function(mm, s, M) as.numeric(predict(mm, newx = M, s = s, type = "link"))
    b_nn <- as.numeric(coef(m_nn, s = "lambda.min"))[-1]
    b_un <- as.numeric(coef(m_un, s = "lambda.min"))[-1]
    ## sum-to-one renormalisation of the SAME non-negative direction
    w_simplex <- if (sum(b_nn) > 0) b_nn / sum(b_nn) else rep(1 / K, K)

    S <- list()   # rule -> list(alt = ., null = .)
    S[["ridge.min (shipped, unconstrained)"]] <-
      list(sc(m_un, "lambda.min", Zf),        sc(m_un, "lambda.min", CAL$Z))
    S[["NNridge.min (probit, b>=0)"]] <-
      list(sc(m_nn, "lambda.min", Zf),        sc(m_nn, "lambda.min", CAL$Z))
    S[["NNridge.1se (probit, b>=0)"]] <-
      list(sc(m_nn, "lambda.1se", Zf),        sc(m_nn, "lambda.1se", CAL$Z))
    S[["NNlogit (probit, b>=0, ~unpenalised)"]] <-
      list(sc(g_nn, lam_end, Zf),             sc(g_nn, lam_end, CAL$Z))
    S[["NNridge.min.SIMPLEX (sum w = 1)"]] <-
      list(as.numeric(Zf %*% w_simplex),      as.numeric(CAL$Z %*% w_simplex))

    ## ---------------- (C) Cauchy-coordinate simplex weights ----------------
    ## INNER leave-one-training-family-out CV to pick (gamma, agg)
    inner <- matrix(NA_real_, nrow = nrow(CFG), ncol = length(trf),
                    dimnames = list(NULL, trf))
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
    ## the same inner CV for the crude top-m subset rule
    innerM <- matrix(NA_real_, nrow = length(MS), ncol = length(trf), dimnames = list(NULL, trf))
    for (gi in seq_along(trf)) {
      g_out <- trf[gi]; a_mean <- ahat_of(setdiff(trf, g_out), "mean")
      for (mi in seq_along(MS)) {
        w <- w_topm(a_mean, MS[mi])
        innerM[mi, gi] <- pow_at(as.numeric(TR[[g_out]]$T %*% w), as.numeric(CAL$T %*% w))
      }
    }
    mbest <- order(-round(rowMeans(innerM), 6), -round(apply(innerM, 1, min), 6), -MS)[1]
    w <- w_topm(ahat_of(trf, "mean"), MS[mbest])
    S[["CCT top-m subset (nested-CV m)"]] <- list(as.numeric(Tf %*% w), as.numeric(CAL$T %*% w))
    sel[[length(sel) + 1L]] <- data.frame(n = n, fold = f, criterion = "topm",
                                          gamma = MS[mbest], agg = "mean",
                                          inner_mean = rowMeans(innerM)[mbest],
                                          inner_worst = apply(innerM, 1, min)[mbest])

    inner_mean <- rowMeans(inner); inner_worst <- apply(inner, 1, min)
    ## select by inner-CV MEAN power; ties -> better worst-case -> smaller gamma
    ord  <- order(-round(inner_mean, 6), -round(inner_worst, 6), CFG$gamma)
    best <- ord[1]
    ## the alternative selection criterion (worst-case), reported separately
    ordw <- order(-round(inner_worst, 6), -round(inner_mean, 6), CFG$gamma)
    bestw <- ordw[1]

    for (tag in c("mean", "worst")) {
      ci <- if (tag == "mean") best else bestw
      a  <- ahat_of(trf, CFG$agg[ci])          # refit ahat on all 11 in-fold families
      w  <- w_of(a, CFG$gamma[ci])
      nm <- sprintf("NNCauchy (nested-CV, select=%s)", tag)
      S[[nm]] <- list(as.numeric(Tf %*% w), as.numeric(CAL$T %*% w))
      sel[[length(sel) + 1L]] <- data.frame(n = n, fold = f, criterion = tag,
                                            gamma = CFG$gamma[ci], agg = CFG$agg[ci],
                                            inner_mean = inner_mean[ci],
                                            inner_worst = inner_worst[ci])
      if (tag == "mean")
        wrows[[length(wrows) + 1L]] <- data.frame(n = n, fold = f,
          feature = make.names(FEATS11), w = as.numeric(w),
          ahat_mean = as.numeric(ahat_of(trf, "mean")),
          b_NNridge = b_nn, b_ridge_unconstrained = b_un)
    }

    ## ---- diagnostic fixed-gamma curve (NOT used for selection) ------------
    a_mean_full <- ahat_of(trf, "mean")
    for (gm in GAMMAS) {
      w <- w_of(a_mean_full, gm)
      S[[sprintf("  diag: NNCauchy gamma=%s agg=mean", format(gm))]] <-
        list(as.numeric(Tf %*% w), as.numeric(CAL$T %*% w))
    }
    for (m in MS) {
      w <- w_topm(a_mean_full, m)
      S[[sprintf("  diag: CCT top-%02d subset", m)]] <-
        list(as.numeric(Tf %*% w), as.numeric(CAL$T %*% w))
    }

    ## ---------------- baselines, same rows / same calibration --------------
    S[["CCT-full (equal weight) BASELINE"]] <-
      list(rowMeans(Tf), rowMeans(CAL$T))
    S[["Zmean (equal weight, probit)"]] <-
      list(rowMeans(Zf), rowMeans(CAL$Z))
    S[["Zmax (max, probit)"]] <-
      list(apply(Zf, 1, max), apply(CAL$Z, 1, max))

    ## per-rep rejection indicators, so candidate-vs-baseline can be compared
    ## PAIRED on the identical simulated datasets (much tighter than unpaired MC)
    rej  <- lapply(S, function(x)
      as.integer(x[[1]] > stats::quantile(x[[2]], 1 - ALPHA, na.rm = TRUE, names = FALSE)))
    base <- rej[["CCT-full (equal weight) BASELINE"]]
    for (r in names(S))
      rows[[length(rows) + 1L]] <- data.frame(
        n = n, rule = r, family = f, power = mean(rej[[r]]),
        d_paired = mean(rej[[r]] - base),
        se_paired = stats::sd(rej[[r]] - base) / sqrt(length(base)))
  }

  POW <- do.call(rbind, rows)
  SEL <- do.call(rbind, sel)
  WT  <- do.call(rbind, wrows)

  ## ---- CERTIFICATE 1: gamma = 0 IS CCT-full ---------------------------------
  d0 <- merge(POW[POW$rule == "  diag: NNCauchy gamma=0 agg=mean", c("family", "power")],
              POW[POW$rule == "CCT-full (equal weight) BASELINE", c("family", "power")],
              by = "family")
  say("\nCERTIFICATE 1 (class contains the baseline): max |power(NNCauchy, gamma=0) - power(CCT-full)| = %.3e",
      max(abs(d0$power.x - d0$power.y)))
  wtest <- rep(1 / K, K)
  say("              max |T_w(w=1/K) - CCT statistic| over the conformal null   = %.3e",
      max(abs(as.numeric(CAL$T %*% wtest) - rowMeans(CAL$T))))

  ## ---- CERTIFICATE 2: sum-to-one is vacuous for power ----------------------
  d1 <- merge(POW[POW$rule == "NNridge.min (probit, b>=0)", c("family", "power")],
              POW[POW$rule == "NNridge.min.SIMPLEX (sum w = 1)", c("family", "power")],
              by = "family")
  say("CERTIFICATE 2 (sum-to-one is a positive rescaling -> power invariant):")
  say("              max |power(NNridge.min) - power(same beta, renormalised to sum 1)| = %.3e",
      max(abs(d1$power.x - d1$power.y)))

  ## ---- summary table --------------------------------------------------------
  wide <- reshape(POW[, c("rule", "family", "power")],
                  idvar = "rule", timevar = "family", direction = "wide")
  names(wide) <- sub("^power\\.", "", names(wide))
  wide$mean_ALT11  <- round(rowMeans(wide[, ALT11, drop = FALSE]), 4)
  wide$worst_ALT11 <- round(apply(wide[, ALT11, drop = FALSE], 1, min), 4)
  wide[, POOL_OLD12] <- round(wide[, POOL_OLD12], 4)
  wide$role <- ifelse(grepl("^  diag", wide$rule), "diagnostic",
               ifelse(grepl("BASELINE", wide$rule), "baseline",
               ifelse(grepl("^NNCauchy", wide$rule), "CANDIDATE",
               ifelse(grepl("^NNridge|^NNlogit", wide$rule), "candidate", "reference"))))
  wide$n <- n
  ## paired-vs-baseline summary columns (same reps, same calibration)
  pag <- do.call(rbind, lapply(split(POW[POW$family %in% ALT11, ],
                                     POW$rule[POW$family %in% ALT11]), function(g)
    data.frame(rule = g$rule[1],
               d_vs_CCTfull = round(mean(g$d_paired), 4),
               se_paired    = round(sqrt(sum(g$se_paired^2)) / nrow(g), 4),
               wins_of_11   = sum(g$d_paired > 0))))
  pag$t_paired <- round(ifelse(pag$se_paired > 0, pag$d_vs_CCTfull / pag$se_paired, 0), 2)
  wide <- merge(wide, pag, by = "rule")
  wide <- wide[order(-wide$mean_ALT11), ]
  wide <- wide[, c("n", "role", "rule", POOL_OLD12, "mean_ALT11", "worst_ALT11",
                   "d_vs_CCTfull", "se_paired", "t_paired", "wins_of_11")]
  utils::write.csv(wide, file.path(OUT, sprintf("cand_nonneg_trainCV_n%d.csv", n)),
                   row.names = FALSE)
  utils::write.csv(WT,  file.path(OUT, sprintf("cand_nonneg_weights_n%d.csv", n)), row.names = FALSE)
  utils::write.csv(SEL, file.path(OUT, sprintf("cand_nonneg_selection_n%d.csv", n)), row.names = FALSE)

  say("\n-- TRAINING-CV (LOTFO) power, probit excluded from the summary --")
  for (i in seq_len(nrow(wide)))
    say("   %-10s %-38s mean %.4f   worst %.4f",
        wide$role[i], wide$rule[i], wide$mean_ALT11[i], wide$worst_ALT11[i])

  bl <- wide[wide$rule == "CCT-full (equal weight) BASELINE", ]
  say("\n-- vs the CCT-full training baseline (mean %.4f / worst %.4f) --",
      bl$mean_ALT11, bl$worst_ALT11)
  for (r in c("NNridge.min (probit, b>=0)", "NNCauchy (nested-CV, select=mean)",
              "NNCauchy (nested-CV, select=worst)", "ridge.min (shipped, unconstrained)")) {
    i <- which(wide$rule == r)
    say("   %-38s d(mean) %+0.4f   d(worst) %+0.4f", r,
        wide$mean_ALT11[i] - bl$mean_ALT11, wide$worst_ALT11[i] - bl$worst_ALT11)
  }

  ## per-family sign test + PAIRED Monte-Carlo error vs the baseline
  say("\n-- PAIRED comparison with CCT-full on the identical simulated datasets --")
  say("   (families are independent draws, so se(mean d) = sqrt(sum se_f^2)/11;")
  say("    the shared conformal-null threshold is common to both rules and ignored)")
  for (r in c("NNridge.min (probit, b>=0)", "NNCauchy (nested-CV, select=mean)",
              "CCT top-m subset (nested-CV m)", "ridge.min (shipped, unconstrained)")) {
    g  <- POW[POW$rule == r & POW$family %in% ALT11, ]
    d  <- g$d_paired
    se <- sqrt(sum(g$se_paired^2)) / length(d)
    say("   %-38s mean d %+0.4f  se %0.4f  (t = %+0.2f)  wins %d/%d families",
        r, mean(d), se, mean(d) / se, sum(d > 0), length(d))
  }

  ## ---- selected hyper-parameters + weight structure -------------------------
  s1 <- SEL[SEL$criterion == "mean", ]
  say("\n-- nested-CV selected tilt (criterion = mean) over the %d outer folds --", nrow(s1))
  say("   gamma chosen: %s", paste(sprintf("%s(%s)", format(s1$gamma), s1$agg), collapse = " "))
  say("   gamma = 0 (i.e. the procedure fell back to EXACT CCT-full) in %d/%d folds",
      sum(s1$gamma == 0), nrow(s1))
  s2 <- SEL[SEL$criterion == "worst", ]
  say("   [criterion = worst] gamma chosen: %s",
      paste(sprintf("%s(%s)", format(s2$gamma), s2$agg), collapse = " "))
  s3 <- SEL[SEL$criterion == "topm", ]
  say("   [crude top-m subset rule] m chosen: %s", paste(format(s3$gamma), collapse = " "))

  ## ---- how much did the non-negativity constraint actually bind? -----------
  nneg <- tapply(WT$b_ridge_unconstrained, WT$feature, function(v) sum(v < 0))
  mag  <- tapply(WT$b_ridge_unconstrained, WT$feature, function(v) mean(v))
  say("\n-- how much does beta >= 0 bind? (unconstrained ridge coefficients, %d folds) --", 12L)
  neg <- names(nneg)[nneg > 0]
  if (!length(neg)) say("   no coordinate is ever negative: the constraint is INACTIVE and NNridge == ridge")
  else for (v in neg)
    say("   %-16s negative in %2d/12 folds, mean coefficient %+0.4f", v, nneg[v], mag[v])
  say("   mean Cauchy weight w_j (nested-CV selected, criterion = mean):")
  ws <- tapply(WT$w, WT$feature, mean); ah <- tapply(WT$ahat_mean, WT$feature, mean)
  for (v in names(sort(-ws)))
    say("     %-16s w %.4f  (equal weight = %.4f)   marginal power ahat = %.3f", v, ws[v], 1 / K, ah[v])
}
say("\ncand_nonneg.R: done")
