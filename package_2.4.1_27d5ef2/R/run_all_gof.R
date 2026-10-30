#' Run a Battery of Goodness-of-Fit Tests at Once
#'
#' @description
#' Runs several goodness-of-fit tests for a binary logistic regression in one
#' call and returns one tidy \code{data.frame}, one row per test. Pass a fitted
#' \code{glm} to run the whole battery; pass \code{(y, predicted_probs)} to run
#' the tests that need only predictions. Each test is wrapped so that a failure of
#' one test never aborts the whole run.
#'
#' @details
#' The currently bundled tests are: \code{Pearson}, \code{Deviance},
#' \code{Osius-Rojek}, \code{McCullagh}, \code{Copas-RSS}, and
#' \code{Information-Matrix} (the White/Orme test) (global / standardized);
#' \code{McCullagh} standardizes the Pearson statistic by its exact conditional
#' moments (Kuss 2002 algorithm);
#' \code{HL} (Hosmer-Lemeshow deciles), \code{HL-equalwidth},
#' \code{Pigeon-Heyse}, and \code{F-test} (the modified Hosmer-Lemeshow
#' F-test: deviance residuals ANOVA-F-tested across deciles) (partition); \code{EF} and \code{EF-normal} (the omnibus
#' Ebrahim-Farrington test with the chi-square and normal references; the normal
#' form reproduces the thesis simulation); \code{DEF.poly2/poly3/stukel}
#' and \code{Stukel} (directed); \code{Tsiatis}, \code{Xie}, and
#' \code{Pulkstenis-Robinson} (covariate-space); the two ensemble rows
#' (\code{Ensemble.Vote(3DEF)} and \code{Ensemble.Univ(3DEF+EF)}) from the Cauchy
#' combination test; and, when \code{include_slow = TRUE}, the opt-in slow tests:
#' \code{le-Cessie}-van Houwelingen smoothing, the GAM-based tests \code{HL-GAM},
#' \code{PR-GAM}, \code{Xie-GAM} (need \pkg{mgcv}; fit an overfit GAM for grouping),
#' \code{Stute-Zhu} (a cumulative-residual parametric-bootstrap test; set the
#' number of reps with \code{control = list("Stute-Zhu" = list(B = ...))}),
#' \code{eHL} (the e-value Hosmer-Lemeshow test, reported as p = min(1, 1/e)), and
#' \code{BAGofT} (the binary-adaptive GOF test; needs the \pkg{BAGofT} package,
#' \code{control = list(BAGofT = list(nsim = ...))}), and \code{Lai-Liu-HL} (Lai
#' & Liu's standardized-power procedure for the Hosmer-Lemeshow test, which has no
#' p-value: it reports the standardized power as the statistic and a randomized
#' accept/reject decision in the \code{Note}; target size via
#' \code{control = list("Lai-Liu-HL" = list(n0 = ..., k = ...))}), and
#' \code{GiViTI} and \code{GiViTI-external} (the GiViTI polynomial calibration
#' test with the internal and external development assumptions; wraps
#' \pkg{givitiR}, run in an isolated \pkg{callr} subprocess so a failure in its
#' compiled dependencies returns \code{NA} rather than aborting the session;
#' set \code{control = list(GiViTI = list(devel = "internal"))}).
#'
#' Notes: \code{Tsiatis} and \code{Xie} cluster the covariate space with k-means
#' (a fixed internal seed, so results are reproducible and the caller's RNG is
#' left untouched). \code{Xie} uses the corrected degrees of freedom
#' \eqn{G - k/2 - 1} with \eqn{k} the number of predictors. \code{Pulkstenis-Robinson}
#' auto-detects the categorical covariate (any factor/character/logical, or a
#' numeric with at most \code{getOption("ebrahim.gof.pr.maxlev", 6)} distinct
#' values); it returns \code{NA} with a note when none is present.
#'
#' Every bundled test reproduces the implementation used in the original thesis
#' simulation: \code{Osius-Rojek} and \code{Stukel} follow \pkg{LogisticDx}'s
#' \code{gof.glm} (Stukel via \code{statmod::glm.scoretest} when \pkg{statmod} is
#' installed), \code{Copas-RSS} follows \pkg{rms}'s gof residual, \code{HL} follows
#' \code{ResourceSelection::hoslem.test}, and the others match their standalone
#' reference functions; all were checked to agree numerically.
#'
#' @param object A fitted binary logistic \code{\link[stats]{glm}}, or a binary
#'   (0/1) response vector \code{y} (then supply \code{predicted_probs}).
#' @param predicted_probs Numeric predicted probabilities; required when
#'   \code{object} is a \code{y} vector.
#' @param X Optional design matrix; lets the directed (DEF) tests run from the
#'   \code{(y, predicted_probs)} form.
#' @param tests Either \code{"all"} (default) or a character vector of test names
#'   to run (e.g. \code{c("EF","DEF.poly3","HL")}).
#' @param G Integer number of groups passed to the grouping tests (default 10).
#' @param include_slow Logical; when \code{TRUE} (the default) the full battery
#'   runs, including the slow tests: le Cessie-van Houwelingen smoothing
#'   (O(n^2)-O(n^3)), the GAM tests, Stute-Zhu, eHL, BAGofT, and GiViTI. Set
#'   \code{FALSE} for a quick run with the fast tests only. A one-time message
#'   notes this whenever slow tests are included.
#' @param parallel Logical; when \code{TRUE}, the resampling loops of the slow
#'   bootstrap tests (\code{Stute-Zhu} and \code{Lai-Liu-HL}) are run on a local
#'   PSOCK cluster via \code{\link[parallel]{parLapply}} (works on all
#'   platforms, including Windows). All other tests are unaffected. The default
#'   \code{FALSE} keeps every loop sequential, exactly as in previous versions.
#' @param ncores Integer; the number of worker processes used when
#'   \code{parallel = TRUE}. The default \code{NULL} uses
#'   \code{max(1, parallel::detectCores() - 1)}. Values below 2 fall back to
#'   the sequential path.
#' @param calibration_plot Logical; when \code{TRUE} and \code{GiViTI} is among
#'   the tests, also compute and draw the GiViTI calibration belt and store it on
#'   the result (retrievable with \code{plot()}). Default \code{FALSE}.
#' @param install One of \code{"ask"} (default), \code{"no"}, or \code{"yes"},
#'   controlling what happens when a test in the run needs an optional package
#'   that is not installed. In an \emph{interactive} session, \code{"ask"} lists
#'   the missing packages and asks before installing, and \code{"yes"} installs
#'   them without asking; \code{"no"} never installs (the test is just skipped
#'   with a note). In a non-interactive session (scripts, \code{R CMD check})
#'   nothing is ever installed, regardless of this setting. See
#'   \code{\link{gof_install_suggests}}.
#' @param control Optional named list of per-test options. Recognized entries:
#'   \code{"Stute-Zhu" = list(B = ...)} (bootstrap replicates);
#'   \code{GiViTI = list(devel = "internal"/"external")};
#'   \code{"Lai-Liu-HL" = list(n0 = ..., k = ..., alpha = ...)}; and
#'   \code{BAGofT = list(...)} which forwards to the binary adaptive test --
#'   \code{nsim} (resampling iterations; default 100), \code{nsplits}, \code{ne}
#'   (the estimation-split size), and the random-forest partitioner's tuning
#'   \code{Kmax} (maximum number of adaptive partition cells), \code{ntree},
#'   \code{nmin}, \code{mtry}, \code{maxnodes}. Example:
#'   \code{list(BAGofT = list(nsim = 200, Kmax = 8, ntree = 500))}.
#'
#' @return A \code{data.frame} (of class \code{gof_battery}) with columns
#'   \code{Test}, \code{Family}, \code{Statistic}, \code{df}, \code{p_value},
#'   and \code{Note}, one row per test. A dedicated \code{print} method shows the
#'   rows grouped by family with formatted p-values and significance flags; the
#'   underlying columns remain available for programmatic use.
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' x <- runif(n, -3, 3)
#' y <- rbinom(n, 1, 1 / (1 + exp(-(0.6 * x))))
#' fit <- glm(y ~ x, family = binomial())
#'
#' ## quick run: the fast tests only
#' run.all.gof(fit, include_slow = FALSE)
#'
#' ## pick specific tests
#' run.all.gof(fit, tests = c("McCullagh", "Osius-Rojek", "HL"))
#'
#' \donttest{
#' ## the full battery (default include_slow = TRUE); the slow tests may need the
#' ## suggested packages mgcv, BAGofT, givitiR and callr. In an interactive
#' ## session run.all.gof() offers to install any that are missing (install =
#' ## "ask"); see also gof_install_suggests(). Use install = "no" to never ask.
#' run.all.gof(fit, install = "no", control = list("Stute-Zhu" = list(B = 50)))
#'
#' ## draw the GiViTI calibration belt (needs givitiR + callr)
#' res <- run.all.gof(fit, tests = c("McCullagh", "GiViTI"),
#'                    calibration_plot = TRUE)
#' plot(res)   # redraw the stored belt
#'
#' ## run the slow bootstrap loops (Stute-Zhu, Lai-Liu-HL) on a PSOCK cluster
#' set.seed(1)
#' run.all.gof(fit, tests = "Stute-Zhu", parallel = TRUE, ncores = 2,
#'             control = list("Stute-Zhu" = list(B = 50)))
#' }
#'
#' @note \strong{Grouped vs sparse forms.} \code{Pearson}, \code{Deviance} and
#'   \code{McCullagh} are reported in two forms: the default (sparse / one-trial)
#'   form and a \code{"(grouped)"} form computed on the distinct covariate
#'   patterns (each a Binomial\eqn{(m_g, P_g)}). The two are identical when every
#'   covariate pattern is unique (fully sparse data, as in the simulation) and
#'   differ only when patterns repeat (\eqn{m_g > 1}). To avoid clutter, the
#'   \code{"(grouped)"} row is shown \emph{only} when it actually differs from the
#'   sparse form (i.e., when some pattern repeats); on fully sparse data it is a
#'   duplicate and is omitted. \code{Osius-Rojek} is always computed on covariate
#'   patterns, matching its classical (LogisticDx) definition.
#'
#'   \strong{Farrington vs EF.} The \emph{original} Farrington (1996) test is a
#'   grouped (covariate-pattern) test. The Ebrahim-Farrington (\code{EF}) test is
#'   its \emph{sparse-data} counterpart: it does not group by covariate pattern
#'   but forms \code{G} data-dependent bins of the predicted risk, so it applies
#'   directly to fully sparse data. Use \code{EF} for sparse binary data; the
#'   grouped Farrington form is appropriate only when covariate patterns repeat.
#'
#'   \strong{Reproducibility of the parallel path.} With \code{parallel = TRUE}
#'   the cluster's random-number streams are initialized with
#'   \code{\link[parallel]{clusterSetRNGStream}}, seeded deterministically from
#'   the session's current RNG state. Two runs from the same
#'   \code{\link{set.seed}} state (and the same \code{ncores}) therefore give
#'   identical bootstrap p-values. Note that the parallel L'Ecuyer-CMRG streams
#'   necessarily differ from the serial RNG stream, so \code{parallel = TRUE}
#'   results differ (within Monte-Carlo error) from \code{parallel = FALSE}
#'   results at the same seed; this is standard and both are valid. Results also
#'   depend on \code{ncores}, because the replicates are split across workers.
#' @seealso \code{\link{ef.gof}}, \code{\link{def.gof}}, \code{\link{def.ensemble.gof}},
#'   \code{\link{gof_install_suggests}}.
#'
#' @references
#' The aggregated tests are due to their original authors; they are provided
#' here for comparison and credited as follows.
#'
#' Farrington CP (1996). "On Assessing Goodness of Fit of Generalized Linear
#' Models to Sparse Data." \emph{Journal of the Royal Statistical Society B},
#' \strong{58}(2), 349--360. \doi{10.1111/j.2517-6161.1996.tb02086.x}
#'
#' Hosmer DW, Lemeshow S (1980). "Goodness of Fit Tests for the Multiple
#' Logistic Regression Model." \emph{Communications in Statistics -- Theory and
#' Methods}, \strong{9}(10), 1043--1069. \doi{10.1080/03610928008827941}
#'
#' McCullagh P (1985). "On the Asymptotic Distribution of Pearson's Statistic in
#' Linear Exponential Family Models." \emph{International Statistical Review},
#' \strong{53}(1), 61--67. \doi{10.2307/1402880}
#'
#' Osius G, Rojek D (1992). "Normal Goodness-of-Fit Tests for Multinomial Models
#' with Large Degrees of Freedom." \emph{Journal of the American Statistical
#' Association}, \strong{87}(420), 1145--1152.
#' \doi{10.1080/01621459.1992.10476271}
#'
#' le Cessie S, van Houwelingen JC (1991). "A Goodness-of-Fit Test for Binary
#' Regression Models, Based on Smoothing Methods." \emph{Biometrics},
#' \strong{47}(4), 1267--1282. \doi{10.2307/2532385}
#'
#' Stukel TA (1988). "Generalized Logistic Models." \emph{Journal of the
#' American Statistical Association}, \strong{83}(402), 426--431.
#' \doi{10.1080/01621459.1988.10478613}
#'
#' Stute W, Zhu LX (2002). "Model Checks for Generalized Linear Models."
#' \emph{Scandinavian Journal of Statistics}, \strong{29}(3), 535--545.
#' \doi{10.1111/1467-9469.00304}
#'
#' Tsiatis AA (1980). "A Note on a Goodness-of-Fit Test for the Logistic
#' Regression Model." \emph{Biometrika}, \strong{67}(1), 250--251.
#' \doi{10.1093/biomet/67.1.250}
#'
#' Xie XJ, Pendergast J, Clarke W (2008). "Increasing the Power: A Practical
#' Approach to Goodness-of-Fit Test for Logistic Regression Models with
#' Continuous Predictors." \emph{Computational Statistics & Data Analysis},
#' \strong{52}(5), 2703--2713. \doi{10.1016/j.csda.2007.09.027}
#'
#' Pulkstenis E, Robinson TJ (2002). "Two Goodness-of-Fit Tests for Logistic
#' Regression Models with Continuous Covariates." \emph{Statistics in Medicine},
#' \strong{21}(1), 79--93. \doi{10.1002/sim.943}
#'
#' Nattino G, Finazzi S, Bertolini G (2014). "A New Calibration Test and a
#' Reappraisal of the Calibration Belt for the Assessment of Prediction Models
#' Based on Dichotomous Outcomes." \emph{Statistics in Medicine}, \strong{33}(14),
#' 2390--2407. \doi{10.1002/sim.6100}
#'
#' Zhang J, Ding J, Yang Y (2021). "Is a Classification Procedure Good Enough? A
#' Goodness-of-Fit Assessment Tool for Classification Learning." \emph{Journal of
#' the American Statistical Association}. \doi{10.1080/01621459.2021.1979010}
#'
#' @importFrom stats fitted predict model.matrix model.frame coef deviance pchisq binomial glm.fit kmeans median dist anova lm pnorm
#' @importFrom utils capture.output
#' @export
run.all.gof <- function(object, predicted_probs = NULL, X = NULL,
                        tests = "all", G = 10, include_slow = TRUE,
                        parallel = FALSE, ncores = NULL,
                        calibration_plot = FALSE, install = c("ask", "no", "yes"),
                        control = list()) {

  install <- match.arg(install)
  ctx <- .gof_context(object, predicted_probs, X, G = G)
  sel <- if (identical(tests, "all")) names(.GOF_REGISTRY) else intersect(tests, names(.GOF_REGISTRY))
  if (length(sel) == 0) stop("run.all.gof: no known tests selected.")

  # Optional PSOCK cluster for the slow bootstrap loops (Stute-Zhu, Lai-Liu-HL).
  # Reproducible: clusterSetRNGStream is seeded deterministically from the
  # session RNG, so same set.seed() + same ncores => identical results. The
  # sequential default path is completely untouched.
  if (isTRUE(parallel) && isTRUE(include_slow) && ctx$has_model &&
      length(intersect(sel, c("Stute-Zhu", "Lai-Liu-HL"))) > 0) {
    nc <- if (!is.null(ncores)) max(1L, as.integer(ncores)[1]) else {
      dc <- parallel::detectCores()
      if (is.na(dc)) 2L else max(1L, dc - 1L)
    }
    if (nc >= 2L) {
      cl <- tryCatch(parallel::makePSOCKcluster(nc), error = function(e) NULL)
      if (is.null(cl)) {
        message("run.all.gof: could not start a PSOCK cluster; ",
                "running sequentially.")
      } else {
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::clusterCall(cl, function(lp) { .libPaths(lp); NULL }, .libPaths())
        parallel::clusterSetRNGStream(cl, iseed = sample.int(.Machine$integer.max, 1L))
        ctx$cl <- cl
      }
    }
  }
  # If a test in this run needs an optional package that isn't installed, offer
  # to install it (interactive sessions only, with confirmation -- see
  # gof_install_suggests). Non-interactive runs are never touched.
  .gof_maybe_install(sel, include_slow, install)
  if (isTRUE(include_slow) &&
      any(vapply(sel, function(nm) isTRUE(.GOF_REGISTRY[[nm]]$slow), logical(1))))
    message("run.all.gof: running the full battery, including the slow tests ",
            "(le-Cessie, the GAM tests, Stute-Zhu, eHL, BAGofT, GiViTI). ",
            "For a quick run with the fast tests only, set include_slow = FALSE.")

  rows <- list(); skipped_model <- FALSE
  for (nm in sel) {
    e <- .GOF_REGISTRY[[nm]]
    if (isTRUE(e$slow) && !include_slow) next
    if (isTRUE(e$needs_model) && !ctx$has_model && is.null(ctx$X)) {
      skipped_model <- TRUE
      rows[[nm]] <- data.frame(Test = nm, Family = e$family, Statistic = NA_real_,
                               df = NA_real_, p_value = NA_real_,
                               Note = "Not applicable: needs a fitted glm", stringsAsFactors = FALSE)
      next
    }
    res <- tryCatch(e$fn(ctx, control[[nm]]),
                    error = function(err) list(Statistic = NA, df = NA, p_value = NA,
                                               Note = paste("error:", conditionMessage(err))))
    rows[[nm]] <- data.frame(Test = nm, Family = e$family,
                             Statistic = .gof_num(res$Statistic), df = .gof_num(res$df),
                             p_value = .gof_num(res$p_value),
                             Note = if (is.null(res$Note)) "" else res$Note,
                             stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)

  # Drop redundant "(grouped)" rows when covariate patterns do not repeat. On
  # fully sparse data (every observation its own pattern) the grouped and
  # per-observation forms of Pearson / deviance / McCullagh coincide exactly, so
  # the "(grouped)" row is an identical duplicate. Keep it only when the two
  # genuinely differ (i.e., some pattern has m_g > 1).
  .eq_gof <- function(x, y) (is.na(x) && is.na(y)) ||
    (!is.na(x) && !is.na(y) && isTRUE(all.equal(x, y, tolerance = 1e-8)))
  gi <- grep(" \\(grouped\\)$", out$Test)
  drop_rows <- integer(0)
  for (k in gi) {
    base <- sub(" \\(grouped\\)$", "", out$Test[k])
    bi <- which(out$Test == base)
    if (length(bi) == 1L && .eq_gof(out$Statistic[k], out$Statistic[bi]) &&
        .eq_gof(out$p_value[k], out$p_value[bi]))
      drop_rows <- c(drop_rows, k)
  }
  if (length(drop_rows)) out <- out[-drop_rows, , drop = FALSE]

  # ensemble rows (only when a model is available and running the full set)
  if (ctx$has_model && identical(tests, "all")) {
    v3 <- tryCatch(def.ensemble.gof(ctx$model, G = G)$p_value, error = function(e) NA_real_)
    vu <- tryCatch(def.ensemble.gof(ctx$model, add_ef = TRUE, G = G)$p_value, error = function(e) NA_real_)
    out <- rbind(out, data.frame(
      Test = c("Ensemble.Vote(3DEF)", "Ensemble.Univ(3DEF+EF)"), Family = "Ensemble",
      Statistic = NA_real_, df = NA_real_, p_value = c(v3, vu),
      Note = "Cauchy combination of the directed tests",
      stringsAsFactors = FALSE))
  }

  if (skipped_model)
    message("Only prediction-based tests were run. Pass the fitted glm (or X) ",
            "to also run the directed and refit-based tests.")
  rownames(out) <- NULL
  class(out) <- c("gof_battery", "data.frame")

  # Optional GiViTI calibration belt: compute it in an isolated subprocess,
  # store it on the result, and draw it.
  if (isTRUE(calibration_plot) && "GiViTI" %in% sel) {
    belt <- .giviti_belt(ctx, control[["GiViTI"]])
    if (is.null(belt)) {
      message("run.all.gof: GiViTI calibration belt unavailable ",
              "(needs the 'givitiR' and 'callr' packages).")
    } else {
      attr(out, "giviti_belt") <- belt
      try(plot(belt), silent = TRUE)
    }
  }
  out
}

# Compute the GiViTI calibration belt object in an isolated callr subprocess
# (same crash-safety as the GiViTI test). Returns NULL if unavailable.
.giviti_belt <- function(ctx, opts = list()) {
  if (!requireNamespace("givitiR", quietly = TRUE) ||
      !requireNamespace("callr", quietly = TRUE))
    return(NULL)
  devel <- if (is.null(opts$devel)) "internal" else opts$devel
  o <- as.numeric(ctx$y); e <- as.numeric(ctx$ph)
  tryCatch(
    callr::r(function(o, e, devel) {
      givitiR::givitiCalibrationBelt(o = o, e = e, devel = devel)
    }, args = list(o = o, e = e, devel = devel), timeout = 120),
    error = function(err) NULL)
}

#' Plot the GiViTI calibration belt from a goodness-of-fit battery
#'
#' Draws the GiViTI calibration belt stored on a \code{\link{run.all.gof}} result
#' that was produced with \code{calibration_plot = TRUE}. The belt shows the
#' fitted calibration curve with a confidence region against the 45-degree line.
#'
#' @param x A \code{gof_battery} object from \code{\link{run.all.gof}}.
#' @param ... Passed to the \pkg{givitiR} plot method.
#' @return \code{x}, invisibly.
#' @importFrom graphics plot
#' @exportS3Method plot gof_battery
plot.gof_battery <- function(x, ...) {
  belt <- attr(x, "giviti_belt")
  if (is.null(belt)) {
    message("No GiViTI calibration belt is stored on this result. ",
            "Re-run run.all.gof(..., calibration_plot = TRUE).")
    return(invisible(NULL))
  }
  if (!requireNamespace("givitiR", quietly = TRUE)) {
    message("Install the 'givitiR' package to draw the calibration belt.")
    return(invisible(NULL))
  }
  plot(belt, ...)
  invisible(x)
}

#' Print a goodness-of-fit battery
#'
#' Formats the \code{\link{run.all.gof}} result as a compact, readable table:
#' rows grouped by test family, p-values shown to four decimals (or scientific
#' for very small values, \code{"-"} when not available), and a significance
#' flag. The object is still a plain \code{data.frame} underneath, so all the
#' raw columns remain available for programmatic use.
#'
#' @param x A \code{gof_battery} object returned by \code{\link{run.all.gof}}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @exportS3Method print gof_battery
print.gof_battery <- function(x, ...) {
  d <- as.data.frame(x)
  fam_order <- c("Global", "Standardized", "Partition", "Covariate-space",
                 "Directed", "Smoothing", "GAM", "Bootstrap", "Calibration", "Ensemble")
  ord <- match(d$Family, fam_order); ord[is.na(ord)] <- 99L
  d <- d[order(ord, seq_len(nrow(d))), , drop = FALSE]

  ## distinct notes -> footnote markers [a], [b], ... (keeps the table narrow)
  notes <- ifelse(is.na(d$Note), "", d$Note)
  uniq  <- unique(notes[notes != ""])
  marks <- letters[seq_along(uniq)]; names(marks) <- uniq
  mk    <- ifelse(notes == "", "", paste0(" [", marks[notes], "]"))

  fp <- vapply(d$p_value, function(p)
    if (is.na(p)) "-" else if (p < 1e-4) sprintf("%.1e", p)
    else formatC(p, format = "f", digits = 4), character(1))
  sig <- vapply(d$p_value, function(p)
    if (is.na(p)) "" else if (p < .001) "***" else if (p < .01) "**"
    else if (p < .05) "*" else if (p < .1) "." else "", character(1))
  fnum <- function(v, dig) vapply(v, function(z)
    if (is.na(z)) "" else formatC(z, format = "g", digits = dig), character(1))
  st <- fnum(d$Statistic, 4); df <- fnum(d$df, 3)

  tname <- paste0(d$Test, mk)
  wT <- max(nchar(tname), nchar("Test"))
  wS <- max(nchar(st), nchar("Statistic"))
  wD <- max(nchar(df), 2L)
  wP <- max(nchar(fp), nchar("p-value")); wG <- 3L
  tw <- 1L + wT + 2L + wS + 1L + wD + 2L + wP + 1L + wG   # table width

  padL <- function(s, w) formatC(s, width = -w, flag = " ")  # left-justified
  padR <- function(s, w) formatC(s, width = w)               # right-justified
  row  <- function(t, s, dd, p, g, tag = "")
    cat(" ", padL(t, wT), "  ", padR(s, wS), " ", padR(dd, wD), "  ",
        padR(p, wP), " ", padL(g, wG), tag, "\n", sep = "")

  nrej <- sum(!is.na(d$p_value) & d$p_value < .05)
  cat(sprintf("\nGoodness-of-fit battery: %d tests  (%d reject at 0.05)\n",
              nrow(d), nrej))
  cat(strrep("=", tw), "\n", sep = "")
  row("Test", "Statistic", "df", "p-value", "")
  fam_now <- ""
  for (i in seq_len(nrow(d))) {
    if (!identical(d$Family[i], fam_now)) {
      fam_now <- d$Family[i]
      lab <- paste0("--- ", fam_now, " ")
      cat(" ", lab, strrep("-", max(0L, tw - nchar(lab) - 1L)), "\n", sep = "")
    }
    row(tname[i], st[i], df[i], fp[i], sig[i])
  }
  cat(strrep("-", tw), "\n", sep = "")
  cat(" Signif.:  *** <.001   ** <.01   * <.05   . <.1\n")
  if (length(uniq)) {
    cat(" Notes:\n")
    for (u in uniq) cat("   [", marks[u], "] ", u, "\n", sep = "")
  }
  invisible(x)
}

# Coerce a possibly-NULL/empty scalar to a numeric (NA if missing).
.gof_num <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)[1]

# Build the one context object every test reads from.
.gof_context <- function(object, predicted_probs = NULL, X = NULL, G = 10) {
  if (inherits(object, "glm")) {
    if (object$family$family != "binomial")
      stop("run.all.gof: the model must be a binomial glm.")
    y  <- as.numeric(object$y)
    ph <- pmin(pmax(as.numeric(stats::fitted(object)), 1e-6), 1 - 1e-6)
    list(y = y, ph = ph, X = stats::model.matrix(object),
         data = tryCatch(stats::model.frame(object), error = function(e) NULL),
         model = object, G = G, n = length(y), has_model = TRUE,
         p = length(stats::coef(object)))
  } else {
    if (!is.numeric(object))
      stop("run.all.gof: 'object' must be a fitted binomial glm or a numeric (0/1) y vector.")
    y <- as.numeric(object)
    if (is.null(predicted_probs))
      stop("run.all.gof: supply 'predicted_probs' when 'object' is not a glm.")
    ph <- pmin(pmax(as.numeric(predicted_probs), 1e-6), 1 - 1e-6)
    if (length(y) != length(ph))
      stop("run.all.gof: 'object' (y) and 'predicted_probs' lengths differ.")
    if (!all(y %in% c(0, 1)))
      stop("run.all.gof: 'object' must be a binary (0/1) vector or a glm.")
    list(y = y, ph = ph, X = X, data = NULL, model = NULL, G = G, n = length(y),
         has_model = FALSE, p = if (is.null(X)) NA_integer_ else ncol(X))
  }
}

# Equal-frequency grouping of the predicted probabilities into G groups.
.gof_groups_ef <- function(ph, G) {
  n <- length(ph)
  pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
}

# Hosmer-Lemeshow statistic for a given grouping; drops empty/degenerate groups.
.gof_hl_stat <- function(y, ph, grp) {
  idx <- split(seq_along(y), grp)
  O  <- vapply(idx, function(I) sum(y[I]),  numeric(1))
  E  <- vapply(idx, function(I) sum(ph[I]), numeric(1))
  ng <- vapply(idx, length, numeric(1))
  keep <- ng > 0 & E > 1e-8 & E < ng - 1e-8
  hl <- sum((O[keep] - E[keep])^2 / (E[keep] * (1 - E[keep] / ng[keep])))
  list(stat = hl, df = sum(keep) - 2)
}

# ---- test wrappers (internal): each returns list(Statistic, df, p_value, Note) ----

gof_pearson <- function(ctx, opts = list()) {
  V  <- ctx$ph * (1 - ctx$ph)
  X2 <- sum((ctx$y - ctx$ph)^2 / V)
  df <- if (ctx$has_model) ctx$model$df.residual else ctx$n - 1
  list(Statistic = X2, df = df, p_value = stats::pchisq(X2, df, lower.tail = FALSE),
       Note = if (ctx$has_model) "" else "df = n-1 (no model)")
}

gof_deviance <- function(ctx, opts = list()) {
  if (ctx$has_model) {
    D <- ctx$model$deviance; df <- ctx$model$df.residual
  } else {
    D  <- -2 * sum(ctx$y * log(ctx$ph) + (1 - ctx$y) * log(1 - ctx$ph)); df <- ctx$n - 1
  }
  list(Statistic = D, df = df, p_value = stats::pchisq(D, df, lower.tail = FALSE),
       Note = if (ctx$has_model) "" else "df = n-1 (no model)")
}

# Osius-Rojek normal-approximation test, matching LogisticDx::gof.glm exactly.
# The statistic is defined on the J DISTINCT COVARIATE PATTERNS, not the raw rows:
# z = (PrG - (J - p)) / sqrt(A1 + RSS1), where PrG is the grouped Pearson statistic,
# A1 = 2*(J - sum(1/n_g)) is the excess-variance correction over the pattern sizes
# n_g, and RSS1 is the residual SS of the WLS regression of (1-2P_g)/(n_g P_g(1-P_g))
# on the pattern design matrix with weights n_g P_g(1-P_g). For fully sparse data
# (every pattern unique, n_g = 1) this reduces to A1 = 0 and J = N, so simulation
# results are unchanged; only data with repeated covariate patterns differ.
gof_osius <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y; p <- ncol(X)
  # aggregate observations sharing an identical covariate pattern (row of X)
  key  <- apply(round(X, 10), 1L, paste, collapse = "\r")
  ug   <- !duplicated(key)
  nm   <- key[ug]
  ng   <- as.numeric(tapply(rep(1, length(y)), key, sum)[nm])   # trials per pattern
  ysum <- as.numeric(tapply(y,  key, sum)[nm])                  # events per pattern
  Pg   <- as.numeric(tapply(ph, key, `[`, 1L)[nm])              # fitted prob per pattern
  Xu   <- X[ug, , drop = FALSE]                                 # one row per pattern
  J    <- length(ng)
  Vg   <- ng * Pg * (1 - Pg)
  PrG  <- sum((ysum - ng * Pg)^2 / Vg)                          # grouped Pearson
  A1   <- 2 * (J - sum(1 / ng))                                 # = 0 only if all n_g = 1
  cvar <- (1 - 2 * Pg) / Vg
  fit  <- tryCatch(stats::lm.wfit(Xu, cvar, Vg), error = function(e) NULL)
  if (is.null(fit))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: weighted least-squares step failed"))
  RSS1 <- sum(Vg * fit$residuals^2)
  varz <- A1 + RSS1
  if (!is.finite(varz) || varz <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: non-positive variance"))
  z <- (PrG - (J - p)) / sqrt(varz)
  list(Statistic = z, df = NA_real_, p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
       Note = if (J < length(y)) sprintf("grouped to %d covariate patterns", J) else "")
}

# Copas (1989) unweighted residual-sum-of-squares test. Binary expansion is
# trivial (one trial per observation). Ported from goflogit (7 tests.R).
gof_copas <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y
  V <- ph * (1 - ph)
  copas    <- sum((y - ph)^2)
  meacopas <- sum(V)
  W   <- diag(V)
  c12 <- 1 - 2 * ph
  M   <- tryCatch(diag(length(ph)) - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X),
                  error = function(e) NULL)
  if (is.null(M))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: singular information matrix"))
  varcopas <- as.numeric(t(c12) %*% M %*% W %*% c12)
  if (!is.finite(varcopas) || varcopas <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: non-positive variance"))
  z <- (copas - meacopas) / sqrt(varcopas)
  list(Statistic = z, df = 1, p_value = stats::pchisq(z^2, 1, lower.tail = FALSE), Note = "")
}

# White/Orme information-matrix (IM) test. Explained sum of squares from
# regressing the Pearson residuals on the auxiliary regressors [sqrt(V) X |
# sqrt(V)(1-2p) X^2]; statistic ~ chi-square_{ncol(X)}. From IM (infromation
# Matrix).R (IMtest_fast); matches the thesis simulation.
gof_im <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y
  if (any(ph < 1e-8) || any(ph > 1 - 1e-8))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: fitted probabilities too extreme (near 0 or 1)"))
  r       <- (y - ph) / sqrt(ph * (1 - ph))
  w_sqrt  <- sqrt(ph * (1 - ph))
  W_aux   <- cbind(w_sqrt * X, (w_sqrt * (1 - 2 * ph)) * (X^2))
  xtx_inv <- tryCatch(solve(crossprod(W_aux)), error = function(e) NULL)
  if (is.null(xtx_inv))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: singular auxiliary matrix"))
  Wr <- crossprod(W_aux, r)
  im <- as.numeric(crossprod(Wr, xtx_inv %*% Wr))
  df <- ncol(X)
  list(Statistic = im, df = df, p_value = stats::pchisq(im, df, lower.tail = FALSE), Note = "")
}

gof_hl <- function(ctx, opts = list()) {
  h <- .gof_hl_stat(ctx$y, ctx$ph, .gof_groups_ef(ctx$ph, ctx$G))
  if (h$df < 1) return(list(Statistic = h$stat, df = h$df, p_value = NA_real_,
                            Note = "too few non-empty groups"))
  list(Statistic = h$stat, df = h$df,
       p_value = stats::pchisq(h$stat, h$df, lower.tail = FALSE), Note = "")
}

gof_hlw <- function(ctx, opts = list()) {
  br <- seq(0, 1, length.out = ctx$G + 1); br[1] <- -Inf; br[length(br)] <- Inf
  grp <- cut(ctx$ph, breaks = br, labels = FALSE)
  h <- .gof_hl_stat(ctx$y, ctx$ph, grp)
  if (h$df < 1) return(list(Statistic = h$stat, df = h$df, p_value = NA_real_,
                            Note = "too few non-empty equal-width bins"))
  list(Statistic = h$stat, df = h$df,
       p_value = stats::pchisq(h$stat, h$df, lower.tail = FALSE), Note = "")
}

# Pigeon-Heyse J2: Hosmer-Lemeshow-type statistic with a per-group variance
# correction phi_k and df = g - 1. Quantile (equal-frequency) groups. From
# pigeonheyse.R; needs only the response and fitted probabilities.
gof_ph_test <- function(ctx, opts = list()) {
  ph <- ctx$ph; y <- ctx$y; g <- ctx$G
  br  <- stats::quantile(ph, probs = seq(0, 1, length.out = g + 1))
  grp <- tryCatch(cut(ph, breaks = br, labels = FALSE, include.lowest = TRUE),
                  error = function(e) NULL)
  if (is.null(grp) || length(unique(grp[!is.na(grp)])) < 2)
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = "could not form groups (ties in fitted probabilities)"))
  idx  <- split(seq_along(y), grp)
  Ok   <- vapply(idx, function(I) sum(y[I]),   numeric(1))
  nk   <- vapply(idx, length,                  numeric(1))
  pbar <- vapply(idx, function(I) mean(ph[I]), numeric(1))
  Vk   <- nk * pbar * (1 - pbar)
  phik <- vapply(idx, function(I) sum(ph[I] * (1 - ph[I])), numeric(1)) / Vk
  ok   <- is.finite(Vk) & Vk > 0 & is.finite(phik) & phik > 0
  J2   <- sum(((Ok[ok] - nk[ok] * pbar[ok])^2 / Vk[ok]) / phik[ok])
  df   <- length(idx) - 1
  list(Statistic = J2, df = df, p_value = stats::pchisq(J2, df, lower.tail = FALSE), Note = "")
}

gof_ef <- function(ctx, opts = list()) {
  r <- ef.gof(ctx$y, ctx$ph, G = ctx$G)            # chisq reference (package default)
  list(Statistic = r$Test_Statistic, df = ctx$G - 2, p_value = r$p_value, Note = "")
}

# EF with the normal reference: reproduces the thesis simulation's farrington_test.
gof_ef_normal <- function(ctx, opts = list()) {
  r <- ef.gof(ctx$y, ctx$ph, G = ctx$G, method = "normal")
  list(Statistic = r$Test_Statistic, df = ctx$G - 2, p_value = r$p_value,
       Note = "normal reference (thesis)")
}

gof_def <- function(ctx, opts = list()) {
  if (!ctx$has_model && is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model or X"))
  b <- if (is.null(opts$basis)) "poly3" else opts$basis
  r <- if (ctx$has_model) def.gof(ctx$model, G = ctx$G, basis = b)
       else suppressWarnings(def.gof(ctx$y, ctx$ph, X = ctx$X, G = ctx$G, basis = b))
  list(Statistic = r$Test_Statistic, df = r$df, p_value = r$p_value, Note = "")
}

# Stukel (1988) two-direction link test, "SstBoth". Matches the thesis simulation
# (LogisticDx::gof.glm), which uses the Rao SCORE test (statmod::glm.scoretest) on
# the sign-split squared-logit directions: the two marginal score-z values are
# squared and summed to a chi-square_2 statistic.
gof_stukel <- function(ctx, opts = list()) {
  if (!ctx$has_model)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  ph  <- ctx$ph; y <- ctx$y; X <- ctx$X
  eta <- as.numeric(stats::predict(ctx$model, type = "link"))
  za  <- 0.5 * eta^2 * (ph >= 0.5)                 # Stukel direction, p >= 0.5
  zb  <- -0.5 * eta^2 * (ph < 0.5)                 # Stukel direction, p < 0.5
  if (requireNamespace("statmod", quietly = TRUE)) {
    # exact match to the thesis simulation (LogisticDx::gof.glm uses glm.scoretest)
    Z <- abs(statmod::glm.scoretest(ctx$model, cbind(za, zb)))
    chi <- sum(Z^2)
  } else {
    za_z <- .gof_score_z(za, y, ph, X)
    zb_z <- .gof_score_z(zb, y, ph, X)
    if (!is.finite(za_z) || !is.finite(zb_z))
      return(list(Statistic = NA, df = NA, p_value = NA, Note = "score test undefined"))
    chi <- za_z^2 + zb_z^2
  }
  list(Statistic = chi, df = 2, p_value = stats::pchisq(chi, 2, lower.tail = FALSE), Note = "")
}

# Rao score-test z for adding one column to a fitted binomial glm. Equals
# statmod::glm.scoretest for the binomial family (dispersion = 1).
.gof_score_z <- function(z, y, ph, X) {
  W   <- ph * (1 - ph)
  U   <- sum(z * (y - ph))
  WzX <- crossprod(X, W * z)
  Vv  <- sum(W * z^2) - as.numeric(t(WzX) %*% solve(crossprod(X, W * X)) %*% WzX)
  if (!is.finite(Vv) || Vv <= 0) return(NA_real_)
  U / sqrt(Vv)
}

# Moore-Penrose pseudo-inverse via SVD (matches MASS::ginv; avoids a dependency).
.gof_ginv <- function(M, tol = sqrt(.Machine$double.eps)) {
  s <- svd(M)
  pos <- s$d > max(tol * s$d[1], 0)
  if (!any(pos)) return(matrix(0, ncol(M), nrow(M)))
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

# k-means with a fixed seed (123) that does not disturb the caller's RNG state.
.gof_kmeans <- function(mat, centers, nstart = 1L) {
  has <- exists(".Random.seed", envir = .GlobalEnv)
  old <- if (has) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(123)
  cl <- stats::kmeans(mat, centers = centers, nstart = nstart)$cluster
  if (has) assign(".Random.seed", old, envir = .GlobalEnv)
  cl
}

# Tsiatis (1980) clustering score test: cluster the covariate space, then score-
# test the cluster indicators added to the model. Ported from Tsiatis.R.
gof_tsiatis <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y
  cov_mat <- X[, -1, drop = FALSE]                       # drop intercept column
  if (ncol(cov_mat) < 1)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "no covariates to cluster"))
  cl <- tryCatch(.gof_kmeans(cov_mat, ctx$G), error = function(e) NULL)
  if (is.null(cl))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "clustering failed"))
  Xc <- stats::model.matrix(~ factor(cl))[, -1, drop = FALSE]
  if (ncol(Xc) < 1)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "only one non-empty cluster"))
  W   <- ph * (1 - ph)
  U   <- colSums(Xc * (y - ph))
  V11 <- crossprod(X, W * X)
  V12 <- crossprod(X, W * Xc)
  V22 <- crossprod(Xc, W * Xc)
  V   <- V22 - t(V12) %*% .gof_ginv(V11) %*% V12
  Tstat <- as.numeric(t(U) %*% .gof_ginv(V) %*% U)
  rankV <- sum(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 1e-8)
  list(Statistic = Tstat, df = rankV,
       p_value = stats::pchisq(Tstat, rankV, lower.tail = FALSE), Note = "")
}

# Xie covariate-space grouped chi-square (own group rule, fractional df). From Xie.R.
gof_xie <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y
  k <- ncol(X) - 1
  cov_mat <- X[, -1, drop = FALSE]
  if (ncol(cov_mat) < 1)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "no covariates to cluster"))
  G <- if (k < 5) 10 else k + 5
  cl <- tryCatch(.gof_kmeans(cov_mat, G, nstart = 25L), error = function(e) NULL)
  if (is.null(cl))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "clustering failed"))
  stat <- 0
  for (I in split(seq_along(y), cl)) {
    ng <- length(I); pbar <- mean(ph[I])
    if (pbar > 0 && pbar < 1)
      stat <- stat + (sum(y[I]) - ng * pbar)^2 / (ng * pbar * (1 - pbar))
  }
  df <- G - k / 2 - 1
  if (df <= 0) return(list(Statistic = stat, df = df, p_value = NA, Note = "df <= 0"))
  list(Statistic = stat, df = df, p_value = stats::pchisq(stat, df, lower.tail = FALSE), Note = "")
}

# Pulkstenis-Robinson: covariate patterns from categorical vars, split by median
# fitted prob, chi-square on the 2*M subgroups. Base-R port of PR_test_only.R.
gof_pr <- function(ctx, opts = list()) {
  if (!ctx$has_model || is.null(ctx$data))
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = "needs a glm model (for categorical covariates)"))
  # model.frame stores the response in column 1; the rest are the model terms.
  cand <- names(ctx$data)[-1]
  maxlev <- getOption("ebrahim.gof.pr.maxlev", 6)
  cat_vars <- if (!is.null(opts$cat_var)) opts$cat_var else
    cand[vapply(cand, function(v) {
      col <- ctx$data[[v]]
      is.factor(col) || is.character(col) || is.logical(col) ||
        (is.numeric(col) && length(unique(col)) <= maxlev)
    }, logical(1))]
  if (length(cat_vars) == 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs a categorical covariate"))

  y <- ctx$y; ph <- ctx$ph
  patt <- do.call(paste, c(lapply(cat_vars, function(v) as.character(ctx$data[[v]])), sep = "_"))
  M <- length(unique(patt))
  lev <- character(length(y))
  for (pp in unique(patt)) {
    ix <- which(patt == pp)
    lev[ix] <- ifelse(ph[ix] <= stats::median(ph[ix]), "low", "high")
  }
  idx <- split(seq_along(y), paste(patt, lev, sep = "::"))
  os <- vapply(idx, function(I) sum(y[I] == 1), numeric(1))
  of <- vapply(idx, function(I) sum(y[I] == 0), numeric(1))
  es <- vapply(idx, function(I) sum(ph[I]), numeric(1))
  ef <- vapply(idx, function(I) sum(1 - ph[I]), numeric(1))
  keep <- es > 0 & ef > 0
  chisq <- sum((os[keep] - es[keep])^2 / es[keep] + (of[keep] - ef[keep])^2 / ef[keep])
  df <- 2 * M - length(cat_vars) - 2
  if (df <= 0)
    return(list(Statistic = chisq, df = df, p_value = NA, Note = "df <= 0 (too few patterns)"))
  list(Statistic = chisq, df = df, p_value = stats::pchisq(chisq, df, lower.tail = FALSE),
       Note = paste0("split on categorical: ", paste(cat_vars, collapse = ", ")))
}

# le Cessie-van Houwelingen smoothed-residual GOF test (general, multivariate).
# Reference: le Cessie, S. & van Houwelingen, H.C. (1995), Biometrics 51:600-614.
# Adapted (with attribution) from the USGS 'smwrStats' package leCessie.test(),
# which is a work of the US federal government (public domain). It builds an
# n-by-n smoothing/kernel matrix, so it is O(n^2)-O(n^3): a Tier-2 ('slow') test.
gof_lecessie <- function(ctx, opts = list()) {
  if (!ctx$has_model || is.null(ctx$data))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  fits   <- ctx$ph
  resids <- ctx$y - fits
  N      <- length(resids)
  covs   <- ctx$data[, -1, drop = FALSE]                 # model.frame minus response
  if (ncol(covs) < 1)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "no covariates"))

  # per-covariate squared distances (lower-triangle vectors), summed across covariates
  dlist <- lapply(covs, function(x) {
    if (is.numeric(x)) {
      as.numeric(0.5 * (stats::dist(scale(x)))^2)
    } else {
      xx <- as.numeric(as.factor(x)); nc <- length(unique(xx))
      if (nc <= 1) rep(0, N * (N - 1) / 2)
      else as.numeric((stats::dist(xx, method = "manhattan") != 0) * nc / (nc - 1))
    }
  })
  dist.mat <- matrix(0, N, N)
  dist.mat[lower.tri(dist.mat)] <- sqrt(rowSums(as.data.frame(dlist)))
  dist.mat <- dist.mat + t(dist.mat)
  bandwidth <- if (!is.null(opts$bandwidth)) opts$bandwidth else mean(dist.mat)
  R.raw <- pmax(1 - dist.mat / bandwidth, 0)
  Q.raw <- sum(as.numeric(resids %*% R.raw) * resids)

  X   <- ctx$X
  mu2 <- fits * (1 - fits)
  hat <- (mu2 * X) %*% solve(crossprod(X, mu2 * X)) %*% t(X)   # V X (X'VX)^{-1} X'
  ## Moment reference for Q = r'Rr with r ~ (I-H)(y-mu):  M = (I-H)' R (I-H).
  ## FIX (2026-07-29): this previously computed (I-H) R (I-H).  In WEIGHTED logistic
  ## regression H = VX(X'VX)^{-1}X' is NOT symmetric (it is in OLS, where the
  ## smwrStats original lives), so the transpose matters: the old form inflated the
  ## null size of the reference (empirical 0.063 vs 0.054 exact-moment at n=1000).
  IH    <- diag(N) - hat
  R.cor <- crossprod(IH, R.raw %*% IH)                         # (I-H)' R (I-H)
  E.Q   <- sum(diag(R.cor) * mu2)
  mu4   <- mu2 * (1 - 3 * mu2)
  VarQ1 <- sum(diag(R.cor)^2 * (mu4 - 3 * mu2^2))
  R.tmp <- R.cor * rep(mu2, each = N)
  VarQ2 <- 2 * sum(diag(R.tmp %*% R.tmp))
  VarQ  <- VarQ1 + VarQ2
  if (!is.finite(VarQ) || VarQ <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: non-positive variance"))
  Test <- Q.raw * 2 * E.Q / VarQ
  df   <- 2 * E.Q^2 / VarQ
  list(Statistic = Test, df = df, p_value = stats::pchisq(Test, df, lower.tail = FALSE), Note = "")
}

# ---- Tier-2 (opt-in, slow) tests ----

# dplyr::ntile equivalent (equal-sized groups, larger groups first).
.gof_ntile <- function(x, n) {
  r <- rank(x, ties.method = "first")
  as.integer(floor((n * (r - 1) / length(x)) + 1))
}

# Observed-vs-expected chi-square over groups (expected from the tested model's ph).
.gof_oe_chisq <- function(y, ph, grp, dfree) {
  idx <- split(seq_along(y), grp)
  o1 <- vapply(idx, function(I) sum(y[I] == 1), numeric(1))
  o0 <- vapply(idx, function(I) sum(y[I] == 0), numeric(1))
  e1 <- vapply(idx, function(I) sum(ph[I]),     numeric(1))
  e0 <- vapply(idx, function(I) sum(1 - ph[I]), numeric(1))
  keep  <- e1 > 0 & e0 > 0
  chisq <- sum((o1[keep] - e1[keep])^2 / e1[keep] + (o0[keep] - e0[keep])^2 / e0[keep])
  if (dfree <= 0)
    return(list(Statistic = chisq, df = dfree, p_value = NA, Note = "df <= 0"))
  list(Statistic = chisq, df = dfree,
       p_value = stats::pchisq(chisq, dfree, lower.tail = FALSE), Note = "")
}

# Fit the overfit GAM (smooth continuous + main categorical + smooth pair
# interactions) and return its fitted probabilities (used for grouping). From GAM.R.
.gof_gam_pi <- function(ctx) {
  if (!requireNamespace("mgcv", quietly = TRUE) || is.null(ctx$data)) return(NULL)
  dat   <- ctx$data
  resp  <- names(dat)[1]
  preds <- names(dat)[-1]
  if (length(preds) == 0) return(NULL)
  maxlev <- getOption("ebrahim.gof.pr.maxlev", 6)
  is_cat <- vapply(preds, function(v) {
    col <- dat[[v]]
    is.factor(col) || is.character(col) || is.logical(col) ||
      (is.numeric(col) && length(unique(col)) <= maxlev)
  }, logical(1))
  cont <- preds[!is_cat]; cats <- preds[is_cat]
  terms <- character(0)
  if (length(cont) > 0) terms <- c(terms, paste0("s(", cont, ")"))
  if (length(cats) > 0) terms <- c(terms, cats)
  if (length(cont) > 1)
    for (i in 1:(length(cont) - 1)) for (j in (i + 1):length(cont))
      terms <- c(terms, paste0("s(", cont[i], ",", cont[j], ")"))
  if (length(terms) == 0) return(NULL)
  fml <- stats::as.formula(paste(resp, "~", paste(terms, collapse = "+")))
  g <- tryCatch(suppressWarnings(mgcv::gam(fml, family = stats::binomial(), data = dat)),
                error = function(e) NULL)
  if (is.null(g)) return(NULL)
  list(pi = as.numeric(stats::predict(g, type = "response")),
       cont = cont, cats = cats, k = length(preds))
}

gof_gam_hl <- function(ctx, opts = list()) {
  if (!ctx$has_model) return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  gf <- .gof_gam_pi(ctx)
  if (is.null(gf)) return(list(Statistic = NA, df = NA, p_value = NA, Note = "install 'mgcv' / no covariates"))
  .gof_oe_chisq(ctx$y, ctx$ph, .gof_ntile(gf$pi, 10), 10 - 2)
}

gof_gam_pr <- function(ctx, opts = list()) {
  if (!ctx$has_model) return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  gf <- .gof_gam_pi(ctx)
  if (is.null(gf)) return(list(Statistic = NA, df = NA, p_value = NA, Note = "install 'mgcv' / no covariates"))
  if (length(gf$cats) == 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs a categorical covariate"))
  patt <- do.call(paste, c(lapply(gf$cats, function(v) as.character(ctx$data[[v]])), sep = "_"))
  M <- length(unique(patt)); pig <- gf$pi
  lev <- character(length(ctx$y))
  for (pp in unique(patt)) {
    ix <- which(patt == pp); lev[ix] <- ifelse(pig[ix] <= stats::median(pig[ix]), "low", "high")
  }
  .gof_oe_chisq(ctx$y, ctx$ph, paste(patt, lev, sep = "::"), 2 * M - length(gf$cats) - 2)
}

gof_gam_xie <- function(ctx, opts = list()) {
  if (!ctx$has_model) return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  gf <- .gof_gam_pi(ctx)
  if (is.null(gf)) return(list(Statistic = NA, df = NA, p_value = NA, Note = "install 'mgcv' / no covariates"))
  G  <- if (gf$k < 5) 10 else gf$k + 5
  cm <- as.data.frame(ctx$data[c(gf$cont, gf$cats)])
  for (v in gf$cats) cm[[v]] <- as.numeric(as.factor(cm[[v]]))
  cl <- tryCatch(.gof_kmeans(scale(as.matrix(cm)), G, nstart = 10L), error = function(e) NULL)
  if (is.null(cl)) return(list(Statistic = NA, df = NA, p_value = NA, Note = "clustering failed"))
  .gof_oe_chisq(ctx$y, ctx$ph, cl, round(G - gf$k / 2 - 2))
}

# Stute-Zhu cumulative-residual statistic (residuals ordered by linear predictor).
.gof_tsz_stat <- function(resid, eta) {
  n <- length(resid)
  (1 / n^2) * sum(cumsum(resid[order(eta)])^2)
}

# Stute-Zhu GOF test via parametric (model-based) bootstrap. From Stute-Zhu(Bootstrap).R.
# Sequential by default; when run.all.gof(parallel = TRUE) put a PSOCK cluster on
# ctx$cl, the bootstrap replicates run via parallel::parLapply (reproducible via
# clusterSetRNGStream). opts$B sets the number of bootstrap reps (default 200).
gof_stutezhu <- function(ctx, opts = list()) {
  if (!ctx$has_model) return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  B   <- if (is.null(opts$B)) 200L else as.integer(opts$B)
  y   <- ctx$y; ph <- ctx$ph; X <- ctx$X
  eta <- as.numeric(stats::predict(ctx$model, type = "link"))
  Tobs <- .gof_tsz_stat(y - ph, eta)
  # One bootstrap replicate, self-contained (base-only closure so PSOCK workers
  # need neither this package's namespace nor any serialized surroundings).
  one_rep <- function(i, probs, Xmat) {
    yb <- stats::rbinom(nrow(Xmat), 1, probs)
    fb <- tryCatch(suppressWarnings(stats::glm.fit(Xmat, yb, family = stats::binomial())),
                   error = function(e) NULL)
    if (is.null(fb) || !isTRUE(fb$converged)) return(NA_real_)
    r <- yb - fb$fitted.values
    e <- as.numeric(Xmat %*% fb$coefficients)
    (1 / length(r)^2) * sum(cumsum(r[order(e)])^2)
  }
  environment(one_rep) <- baseenv()
  Tb <- if (!is.null(ctx$cl)) {
    unlist(parallel::parLapply(ctx$cl, seq_len(B), one_rep, probs = ph, Xmat = X))
  } else {
    vapply(seq_len(B), one_rep, numeric(1), probs = ph, Xmat = X)
  }
  list(Statistic = Tobs, df = NA_real_, p_value = mean(Tb >= Tobs, na.rm = TRUE),
       Note = paste0(B, " bootstrap reps",
                     if (!is.null(ctx$cl))
                       paste0(" (parallel, ", length(ctx$cl), " workers)") else ""))
}

# --- eHL: e-value Hosmer-Lemeshow (Henzi et al. 2024) ---
# Adapted (base-R reimplementation, with attribution) from the marius-cp/eHL
# repository. Splits the data, fits isotonic recalibration on the training half,
# and accumulates an e-value on the test half; p = min(1, 1 / mean(e-values)).

.gof_ehl_interp <- function(x_min, x_max, q, xout) {
  x_min[which.min(x_min)] <- 0
  x_max[which.max(x_max)] <- 1
  p  <- c(x_max, x_min)
  qq <- c(q, q)
  o  <- order(p)
  stats::approx(x = p[o], y = qq[o], xout = xout, method = "linear", ties = "ordered")$y
}

.gof_ehl <- function(y, P, boot = 10L, s = 0.5) {
  n  <- length(y)
  ev <- numeric(boot)
  for (t in seq_len(boot)) {
    idx <- sample(n, floor(n * s), replace = FALSE)
    o1  <- order(P[idx]);  P1 <- P[idx][o1];  y1 <- y[idx][o1]      # training, sorted
    o2  <- order(P[-idx]); P2 <- P[-idx][o2]; y2 <- y[-idx][o2]     # test, sorted
    iso <- stats::isoreg(P1, y1)
    iK  <- iso$iKnots
    red <- iK[!duplicated(iso$yf[iK], fromLast = TRUE)]
    bin <- rep.int(seq_along(red), times = diff(c(0, red)))
    sp  <- split(seq_along(P1), bin)
    nb  <- length(sp)
    x_min <- x_max <- q <- numeric(nb)
    for (b in seq_len(nb)) {
      I  <- sp[[b]]; sz <- length(I); sy <- sum(y1[I])
      x_min[b] <- min(P1[I]); x_max[b] <- max(P1[I])
      hatpi    <- sy / sz
      q[b]     <- if (hatpi %in% c(0, 1)) (sy + 0.5) / (sz + 1) else hatpi
    }
    qte   <- .gof_ehl_interp(x_min, x_max, q, P2)
    E     <- (qte^y2 * (1 - qte)^(1 - y2)) / (P2^y2 * (1 - P2)^(1 - y2))
    ev[t] <- prod(E)
  }
  mean(ev)
}

gof_ehl <- function(ctx, opts = list()) {
  boot <- if (is.null(opts$boot)) 10L else as.integer(opts$boot)
  HLe  <- tryCatch(.gof_ehl(ctx$y, ctx$ph, boot = boot, s = 0.5), error = function(e) NA_real_)
  if (!is.finite(HLe))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: e-value computation failed"))
  list(Statistic = HLe, df = NA_real_, p_value = min(1, 1 / HLe),
       Note = "e-value test (reported as p = min(1, 1/e))")
}

# BAGofT (binary-adaptive GOF test) via the BAGofT package. The random-forest
# partitioner needs at least two predictors; for a single-predictor model we add
# a constant helper column to the data (not the formula), the workaround
# documented in Kuss (2002) / the thesis, so the test runs instead of erroring.
gof_bagoft <- function(ctx, opts = list()) {
  if (!ctx$has_model || is.null(ctx$data))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs a fitted glm"))
  if (!requireNamespace("BAGofT", quietly = TRUE))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: install the 'BAGofT' package"))
  nsim  <- if (is.null(opts$nsim)) 100L else as.integer(opts$nsim)
  dat   <- ctx$data
  attr(dat, "terms") <- NULL                       # a model.frame's terms attr breaks BAGofT
  preds <- names(dat)[-1]
  added_const <- FALSE
  if (length(preds) < 2) {                          # single predictor: add a constant column
    dat[[".bagoft_const"]] <- 1
    preds <- c(preds, ".bagoft_const")
    added_const <- TRUE
  }
  link  <- ctx$model$family$link
  # Adaptive random-forest partitioner: pass through any tuning the caller sets
  # via control = list(BAGofT = list(...)). parRF knobs: Kmax (max number of
  # partition cells - the granularity of the adaptive split), ntree, nmin, mtry,
  # maxnodes. BAGofT() knobs: nsim (resampling iterations), nsplits, ne (size of
  # the estimation split).
  rf_args <- list(parVar = preds)
  for (a in c("Kmax", "nmin", "ntree", "mtry", "maxnodes"))
    if (!is.null(opts[[a]])) rf_args[[a]] <- opts[[a]]
  parFun <- do.call(BAGofT::parRF, rf_args)
  bag_args <- list(testModel = BAGofT::testGlmBi(formula = stats::formula(ctx$model), link = link),
                   parFun = parFun, data = dat, nsim = nsim)
  if (!is.null(opts$nsplits)) bag_args$nsplits <- as.integer(opts$nsplits)
  if (!is.null(opts$ne))      bag_args$ne      <- opts$ne
  r <- NULL                                          # silence BAGofT's per-sim console chatter
  invisible(suppressMessages(utils::capture.output(
    r <- tryCatch(suppressWarnings(do.call(BAGofT::BAGofT, bag_args)),
                  error = function(e) NULL))))
  if (is.null(r) || is.null(r$p.value))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: BAGofT computation failed"))
  extra <- paste0(if (!is.null(opts$Kmax)) paste0(", Kmax=", opts$Kmax) else "",
                  if (!is.null(opts$nsplits)) paste0(", nsplits=", opts$nsplits) else "")
  list(Statistic = NA_real_, df = NA_real_, p_value = as.numeric(r$p.value),
       Note = paste0("adaptive RF partition; nsim=", nsim, extra,
                     if (added_const) "; constant column added (single predictor)" else ""))
}

# McCullagh (1985) exact-conditional-moments standardization of the Pearson
# statistic (SAS GOFLOGIT / Kuss 2002 algorithm, ungrouped binary case). The
# Pearson X^2 is standardized by its conditional mean and variance given the
# fitted coefficients, then Z is referred to the normal. Verified to reproduce
# the thesis low-birth-weight result (p = 0.937) to machine precision.
gof_mccullagh <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y; n <- length(y); p <- ncol(X)
  V  <- ph * (1 - ph)
  X2 <- sum((y - ph)^2 / V)
  W  <- diag(V)
  inv <- tryCatch(solve(t(X) %*% W %*% X), error = function(e) NULL)
  if (is.null(inv))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: singular information matrix"))
  H    <- X %*% inv %*% t(X)               # X (X'WX)^{-1} X'
  h    <- diag(H)                          # leverages
  u    <- (1 - 2 * ph) / V
  uhat <- as.numeric(H %*% W %*% u)
  k2   <- V; k3 <- V * (1 - 2 * ph)
  k4   <- ph - 7 * ph^2 + 12 * ph^3 - 6 * ph^4
  E    <- (n - p) - 0.5 * sum(k4 * h / k2) + 0.5 * sum(uhat * k3 * h)
  RSSu <- as.numeric(t(u) %*% (W - W %*% H %*% W) %*% u)
  Var  <- (1 - p / n) * RSSu               # ungrouped: 2*sum((m-1)/m) = 0
  if (!is.finite(Var) || Var <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: non-positive variance"))
  Z <- (X2 - E) / sqrt(Var)
  list(Statistic = Z, df = NA_real_, p_value = stats::pnorm(Z, lower.tail = FALSE), Note = "")
}

# ---------------------------------------------------------------------------
# GROUPED (covariate-pattern) variants of the global / standardized tests.
#
# Pearson, deviance and McCullagh are classically defined on the J distinct
# COVARIATE PATTERNS, each a Binomial(m_g, P_g). The default rows above use the
# sparse / one-trial form (m_g = 1), which is what the sparse-data simulation
# needs and is identical to the grouped form when every pattern is unique. On
# data with REPEATED covariate patterns (m_g > 1) the two differ, so these
# "(grouped)" variants are reported alongside the ungrouped ones. For fully
# sparse data they coincide with the rows above.
#
# NOTE on Farrington vs EF: the *original* Farrington (1996) test is a GROUPED
# (covariate-pattern, m_g > 1) test. The Ebrahim-Farrington (EF) test in this
# package is deliberately its SPARSE-data counterpart: it does NOT group by
# covariate pattern but forms G data-dependent bins of the predicted risk, so it
# applies directly to fully sparse (one-observation-per-pattern) data. EF and the
# grouped Farrington therefore answer the same question by different partitions.
# ---------------------------------------------------------------------------
.gof_patterns <- function(X, ph, y) {
  key <- apply(round(X, 10), 1L, paste, collapse = "\r")
  ug  <- !duplicated(key); nm <- key[ug]
  list(m  = as.numeric(tapply(rep(1, length(y)), key, sum)[nm]),
       yg = as.numeric(tapply(y,  key, sum)[nm]),
       Pg = as.numeric(tapply(ph, key, `[`, 1L)[nm]),
       Xu = X[ug, , drop = FALSE], J = sum(ug))
}

gof_pearson_grouped <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X"))
  g <- .gof_patterns(ctx$X, ctx$ph, ctx$y); p <- ncol(ctx$X)
  Vg <- g$m * g$Pg * (1 - g$Pg)
  X2 <- sum((g$yg - g$m * g$Pg)^2 / Vg); df <- g$J - p
  list(Statistic = X2, df = df, p_value = stats::pchisq(X2, df, lower.tail = FALSE),
       Note = sprintf("grouped to %d covariate patterns", g$J))
}

gof_deviance_grouped <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X"))
  g <- .gof_patterns(ctx$X, ctx$ph, ctx$y); p <- ncol(ctx$X)
  ys <- pmax(g$yg, 1e-10); ms <- pmax(g$m - g$yg, 1e-10)
  D  <- 2 * sum(g$yg * log(ys / (g$m * g$Pg)) + (g$m - g$yg) * log(ms / (g$m * (1 - g$Pg))))
  df <- g$J - p
  list(Statistic = D, df = df, p_value = stats::pchisq(D, df, lower.tail = FALSE),
       Note = sprintf("grouped to %d covariate patterns", g$J))
}

gof_mccullagh_grouped <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X"))
  g <- .gof_patterns(ctx$X, ctx$ph, ctx$y); p <- ncol(ctx$X)
  m <- g$m; P <- g$Pg; yg <- g$yg; Xu <- g$Xu; J <- g$J; Ntot <- sum(m)
  V  <- m * P * (1 - P)
  X2 <- sum((yg - m * P)^2 / V)
  W  <- diag(V)
  inv <- tryCatch(solve(t(Xu) %*% W %*% Xu), error = function(e) NULL)
  if (is.null(inv))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: singular information matrix"))
  H <- Xu %*% inv %*% t(Xu); h <- diag(H)
  u <- (1 - 2 * P) / V; uhat <- as.numeric(H %*% W %*% u)
  k2 <- V; k3 <- V * (1 - 2 * P); k4 <- m * P * (1 - P) * (1 - 6 * P * (1 - P))
  E    <- (J - p) - 0.5 * sum(k4 * h / k2) + 0.5 * sum(uhat * k3 * h)
  RSSu <- as.numeric(t(u) %*% (W - W %*% H %*% W) %*% u)
  Var  <- 2 * sum((m - 1) / m) + (1 - p / Ntot) * RSSu   # grouped: 2*sum((m-1)/m) is nonzero
  if (!is.finite(Var) || Var <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: non-positive variance"))
  Z <- (X2 - E) / sqrt(Var)
  list(Statistic = Z, df = NA_real_, p_value = stats::pnorm(Z, lower.tail = FALSE),
       Note = sprintf("grouped to %d covariate patterns", J))
}

# GiViTI calibration test (Nattino, Finazzi & Bertolini): forward-selects a
# polynomial calibration model and tests it against the identity, using a
# selection-aware null distribution. Wraps givitiR::givitiCalibrationTest, run in
# an isolated callr subprocess so a crash in givitiR's compiled dependencies
# (alabama, rootSolve) returns NA instead of aborting the session. devel defaults
# to "internal" (the model was fit on these data). Verified against the thesis
# low-birth-weight result (internal p = 0.586).
gof_giviti <- function(ctx, opts = list()) {
  if (!requireNamespace("givitiR", quietly = TRUE))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: install the 'givitiR' package"))
  if (!requireNamespace("callr", quietly = TRUE))
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = "Not run: install the 'callr' package (isolates givitiR)"))
  devel <- if (is.null(opts$devel)) "internal" else opts$devel
  o <- as.numeric(ctx$y); e <- as.numeric(ctx$ph)
  pv <- tryCatch(
    callr::r(function(o, e, devel) {
      givitiR::givitiCalibrationTest(o = o, e = e, devel = devel)$p.value
    }, args = list(o = o, e = e, devel = devel), timeout = 120),
    error = function(err) NA_real_)
  if (length(pv) != 1 || !is.finite(pv))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: givitiR calibration test failed"))
  list(Statistic = NA_real_, df = NA_real_, p_value = pv,
       Note = paste0("calibration belt; devel=", devel))
}

# Lai & Liu (2018) standardized-power procedure for the Hosmer-Lemeshow test.
# Reference: Lai, X. & Liu, L. (2018), "A simple test procedure in standardizing
# the power of Hosmer-Lemeshow test in large data sets", J. Statist. Comput.
# Simul. It resamples (stratified by outcome, with replacement) to a target size
# n0, refits the model, computes the HL statistic each time, and estimates the
# rejection rate ("standardized power") at n0. The final decision is randomized:
# reject H0 if U(0,1) < power. There is no p-value: the statistic is the power
# and the accept/reject decision is reported in the Note. From Hosmer Bootstrap.R.
gof_lailiu <- function(ctx, opts = list()) {
  if (!ctx$has_model || is.null(ctx$data))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  if (!requireNamespace("ResourceSelection", quietly = TRUE))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "install 'ResourceSelection'"))
  dat <- ctx$data; attr(dat, "terms") <- NULL
  y   <- dat[[1]]
  n   <- length(y)
  n0    <- if (is.null(opts$n0))    n    else as.integer(opts$n0)   # target sample size
  k     <- if (is.null(opts$k))     200L else as.integer(opts$k)    # resamples
  alpha <- if (is.null(opts$alpha)) 0.05 else opts$alpha
  G     <- ctx$G
  fml   <- stats::formula(ctx$model)
  ev_grp <- dat[y == 1, , drop = FALSE]; ne_grp <- dat[y == 0, , drop = FALSE]
  m <- nrow(ev_grp); nn <- nrow(ne_grp)
  if (m < 2 || nn < 2)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "too few events/non-events"))
  m0  <- max(1, min(round(n0 * (m / n)), m))
  n0e <- max(1, min(n0 - m0, nn))
  crit <- stats::qchisq(1 - alpha, df = G - 2)
  # One resample, self-contained (base + ResourceSelection only) so it can run
  # either sequentially or on the PSOCK cluster run.all.gof(parallel = TRUE)
  # puts on ctx$cl (reproducible via clusterSetRNGStream).
  one_rep <- function(i, ev_grp, ne_grp, m, nn, m0, n0e, fml, G) {
    rs <- rbind(ev_grp[sample(seq_len(m),  m0,  replace = TRUE), , drop = FALSE],
                ne_grp[sample(seq_len(nn), n0e, replace = TRUE), , drop = FALSE])
    mb <- tryCatch(suppressWarnings(stats::glm(fml, data = rs, family = stats::binomial())),
                   error = function(e) NULL)
    if (is.null(mb) || !isTRUE(mb$converged)) return(NA_real_)
    tryCatch(suppressWarnings(as.numeric(
      ResourceSelection::hoslem.test(mb$y, stats::fitted(mb), g = G)$statistic)),
      error = function(e) NA_real_)
  }
  environment(one_rep) <- baseenv()
  hl <- if (!is.null(ctx$cl)) {
    unlist(parallel::parLapply(ctx$cl, seq_len(k), one_rep, ev_grp = ev_grp,
                               ne_grp = ne_grp, m = m, nn = nn, m0 = m0,
                               n0e = n0e, fml = fml, G = G))
  } else {
    vapply(seq_len(k), one_rep, numeric(1), ev_grp = ev_grp, ne_grp = ne_grp,
           m = m, nn = nn, m0 = m0, n0e = n0e, fml = fml, G = G)
  }
  valid <- hl[!is.na(hl)]
  if (length(valid) == 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "no valid resamples"))
  power    <- mean(valid >= crit)
  decision <- if (stats::runif(1) < power) "REJECT H0 (lack of fit)" else "fail to reject H0"
  list(Statistic = power, df = NA_real_, p_value = NA_real_,
       Note = paste0("standardized power=", round(power, 3),
                     " (n0=", n0, "); decision: ", decision))
}

# Hosmer-Lemeshow F-test (the "modified H&L" of LogisticDx::gof.glm): sort by
# fitted probability, cut into G equal-frequency groups, and ANOVA-F-test the
# deviance residuals against the group factor. Matches the thesis F_test row.
gof_ftest <- function(ctx, opts = list()) {
  y <- ctx$y; ph <- ctx$ph
  grp <- factor(.gof_groups_ef(ph, ctx$G))
  if (nlevels(grp) < 2)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: fewer than 2 groups"))
  eps <- 1e-10
  dr  <- sign(y - ph) * sqrt(pmax(0, -2 * (y * log(pmax(ph, eps)) +
                                           (1 - y) * log(pmax(1 - ph, eps)))))
  av <- tryCatch(stats::anova(stats::lm(dr ~ grp)), error = function(e) NULL)
  if (is.null(av) || !("grp" %in% rownames(av)))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: ANOVA failed"))
  list(Statistic = av["grp", "F value"], df = av["grp", "Df"],
       p_value = av["grp", "Pr(>F)"], Note = "deviance residuals ~ groups (ANOVA F)")
}

# Registry of the bundled tests (test wrappers are internal, not exported).
.GOF_REGISTRY <- list(
  "Pearson"             = list(fn = gof_pearson,          family = "Global",       needs_model = FALSE, slow = FALSE),
  "Pearson (grouped)"   = list(fn = gof_pearson_grouped,  family = "Global",       needs_model = TRUE,  slow = FALSE),
  "Deviance"            = list(fn = gof_deviance,         family = "Global",       needs_model = FALSE, slow = FALSE),
  "Deviance (grouped)"  = list(fn = gof_deviance_grouped, family = "Global",       needs_model = TRUE,  slow = FALSE),
  "Osius-Rojek"         = list(fn = gof_osius,            family = "Standardized", needs_model = TRUE,  slow = FALSE),
  "McCullagh"           = list(fn = gof_mccullagh,        family = "Standardized", needs_model = TRUE,  slow = FALSE),
  "McCullagh (grouped)" = list(fn = gof_mccullagh_grouped, family = "Standardized", needs_model = TRUE, slow = FALSE),
  "Copas-RSS"           = list(fn = gof_copas,            family = "Standardized", needs_model = TRUE,  slow = FALSE),
  "Information-Matrix" = list(fn = gof_im,  family = "Global",       needs_model = TRUE,  slow = FALSE),
  "HL"            = list(fn = gof_hl,       family = "Partition",    needs_model = FALSE, slow = FALSE),
  "HL-equalwidth" = list(fn = gof_hlw,      family = "Partition",    needs_model = FALSE, slow = FALSE),
  "Pigeon-Heyse"  = list(fn = gof_ph_test,  family = "Partition",    needs_model = FALSE, slow = FALSE),
  "F-test"        = list(fn = gof_ftest,    family = "Partition",    needs_model = FALSE, slow = FALSE),
  "EF"            = list(fn = gof_ef,       family = "Standardized", needs_model = FALSE, slow = FALSE),
  "EF-normal"     = list(fn = gof_ef_normal, family = "Standardized", needs_model = FALSE, slow = FALSE),
  "DEF.poly2"     = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "poly2")),  family = "Directed", needs_model = TRUE, slow = FALSE),
  "DEF.poly3"     = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "poly3")),  family = "Directed", needs_model = TRUE, slow = FALSE),
  "DEF.stukel"    = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "stukel")), family = "Directed", needs_model = TRUE, slow = FALSE),
  "Stukel"        = list(fn = gof_stukel,   family = "Directed",     needs_model = TRUE,  slow = FALSE),
  "Tsiatis"             = list(fn = gof_tsiatis, family = "Covariate-space", needs_model = TRUE, slow = FALSE),
  "Xie"                 = list(fn = gof_xie,     family = "Covariate-space", needs_model = TRUE, slow = FALSE),
  "Pulkstenis-Robinson" = list(fn = gof_pr,      family = "Covariate-space", needs_model = TRUE, slow = FALSE),
  "le-Cessie"           = list(fn = gof_lecessie, family = "Smoothing",      needs_model = TRUE, slow = TRUE),
  "HL-GAM"              = list(fn = gof_gam_hl,    family = "GAM",            needs_model = TRUE, slow = TRUE),
  "PR-GAM"              = list(fn = gof_gam_pr,    family = "GAM",            needs_model = TRUE, slow = TRUE),
  "Xie-GAM"             = list(fn = gof_gam_xie,   family = "GAM",            needs_model = TRUE, slow = TRUE),
  "Stute-Zhu"           = list(fn = gof_stutezhu,  family = "Bootstrap",      needs_model = TRUE,  slow = TRUE),
  "eHL"                 = list(fn = gof_ehl,        family = "Calibration",    needs_model = FALSE, slow = TRUE),
  "GiViTI"              = list(fn = gof_giviti,      family = "Calibration",    needs_model = FALSE, slow = TRUE),
  "GiViTI-external"     = list(fn = function(ctx, opts) gof_giviti(ctx, list(devel = "external")),
                                                       family = "Calibration",  needs_model = FALSE, slow = TRUE),
  "BAGofT"              = list(fn = gof_bagoft,     family = "Bootstrap",      needs_model = TRUE,  slow = TRUE),
  "Lai-Liu-HL"          = list(fn = gof_lailiu,     family = "Bootstrap",      needs_model = TRUE,  slow = TRUE)
)
