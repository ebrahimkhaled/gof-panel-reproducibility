#' Covariate-Space Directed Ebrahim-Farrington (CDEF) Goodness-of-Fit Test
#'
#' @description
#' A directed goodness-of-fit test for binary logistic regression whose direction
#' lives in \emph{covariate space} (functions of the predictors) rather than in
#' fitted-probability space like \code{\link{def.gof}}. It projects the
#' standardized residuals onto a covariate-space basis (polynomials and pairwise
#' products, natural splines, or a combination that also includes
#' fitted-probability bends) and calibrates the quadratic form with the
#' Farrington estimation-adjusted projection, exactly as in \code{def.gof}. This
#' makes it sensitive to omitted interactions and to local / oscillatory
#' departures that fitted-probability grouping can miss.
#'
#' @details
#' Let \eqn{\tilde r_i=(y_i-\hat p_i)/\sqrt{\hat p_i(1-\hat p_i)}} be the
#' standardized residuals and \eqn{Z} a covariate-space basis matrix. The
#' statistic is \eqn{S=(Z'\tilde r)'(Z'Z)^{-1}(Z'\tilde r)}, whose null
#' distribution is a weighted sum of \eqn{\chi^2_1} variables with weights the
#' eigenvalues of \eqn{(Z'Z)^{-1}Z'\Omega Z}, where
#' \eqn{\Omega=I-V^{1/2}X(X'VX)^{-1}X'V^{1/2}} adjusts for estimating
#' \eqn{\hat\beta}. The p-value uses a Satterthwaite scaled-\eqn{\chi^2}
#' approximation (default) or Imhof's method (\pkg{CompQuadForm}). Rank-deficient
#' bases are reduced automatically.
#'
#' @param object A fitted binary logistic \code{\link[stats]{glm}}, or a binary
#'   (0/1) response vector \code{y} (then supply \code{predicted_probs} and
#'   \code{X}).
#' @param predicted_probs Numeric predicted probabilities; required when
#'   \code{object} is a \code{y} vector.
#' @param X Design/covariate matrix (with or without an intercept column);
#'   required when \code{object} is a \code{y} vector. Ignored when \code{object}
#'   is a glm.
#' @param basis One of \code{"poly"} (squares, cubes, pairwise products),
#'   \code{"spline"} (natural cubic splines per covariate plus a pairwise term;
#'   needs \pkg{splines}), or \code{"combined"} (covariate polynomials plus
#'   fitted-probability bends).
#' @param method One of \code{"satterthwaite"} (default) or \code{"imhof"}.
#'
#' @return A one-row \code{data.frame} with \code{Test}, \code{Basis},
#'   \code{Test_Statistic}, \code{df}, \code{Method}, and \code{p_value}.
#'
#' @references
#' Farrington, C. P. (1996). On Assessing Goodness of Fit of Generalized Linear
#' Models to Sparse Data. \emph{JRSS-B} 58(2), 349-360.
#'
#' @seealso \code{\link{def.gof}}, \code{\link{ef.gof}}, \code{\link{run.all.gof}}.
#'
#' @examples
#' set.seed(1)
#' n <- 600; x1 <- runif(n, -3, 3); x2 <- rnorm(n)
#' # truth has an omitted interaction; fit the additive model
#' y <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 0.4 * x1 * x2))
#' fit <- glm(y ~ x1 + x2, family = binomial())
#' cdef.gof(fit)                    # covariate-space directed test (poly basis)
#' cdef.gof(fit, basis = "spline")  # for local / oscillatory misfit
#'
#' @importFrom stats fitted model.matrix qlogis sd quantile setNames
#' @export
cdef.gof <- function(object, predicted_probs = NULL, X = NULL,
                     basis  = c("poly", "spline", "combined"),
                     method = c("satterthwaite", "imhof")) {
  basis  <- match.arg(basis)
  method <- match.arg(method)

  if (inherits(object, "glm")) {
    if (object$family$family != "binomial")
      stop("'object' must be a binomial glm (or pass y, predicted_probs, X).")
    y  <- as.numeric(object$y)
    ph <- pmin(pmax(as.numeric(stats::fitted(object)), 1e-6), 1 - 1e-6)
    X  <- stats::model.matrix(object)
  } else {
    if (!is.numeric(object)) stop("'object' must be a binomial glm or a 0/1 y vector.")
    y <- as.numeric(object)
    if (is.null(predicted_probs) || is.null(X))
      stop("Provide both 'predicted_probs' and 'X' when 'object' is a y vector.")
    ph <- pmin(pmax(as.numeric(predicted_probs), 1e-6), 1 - 1e-6)
    X  <- as.matrix(X)
    if (!any(apply(X, 2, function(z) all(z == 1)))) X <- cbind(`(Intercept)` = 1, X)
  }
  if (!all(y %in% c(0, 1))) stop("CDEF needs a binary (0/1) response.")
  V  <- ph * (1 - ph)
  r  <- (y - ph) / sqrt(V)
  Xc <- X[, !apply(X, 2, function(z) all(z == 1)), drop = FALSE]   # drop intercept col(s)
  if (ncol(Xc) < 1) stop("No covariates available to build a covariate-space basis.")

  poly_cov <- function() {
    cols <- list()
    for (j in seq_len(ncol(Xc))) { cols[[length(cols) + 1L]] <- Xc[, j]^2
                                   cols[[length(cols) + 1L]] <- Xc[, j]^3 }
    if (ncol(Xc) >= 2) for (a in 1:(ncol(Xc) - 1)) for (b in (a + 1):ncol(Xc))
      cols[[length(cols) + 1L]] <- Xc[, a] * Xc[, b]
    do.call(cbind, cols)
  }
  phat_shape <- function() { e <- stats::qlogis(ph); cbind(e^2, e^3) }

  if (basis == "poly") {
    Z <- poly_cov()
  } else if (basis == "spline") {
    if (!requireNamespace("splines", quietly = TRUE))
      stop("basis = 'spline' requires the 'splines' package.")
    parts <- lapply(seq_len(ncol(Xc)), function(j) {
      xj <- Xc[, j]; if (length(unique(xj)) > 5) splines::ns(xj, df = 4) else NULL })
    parts <- parts[!vapply(parts, is.null, logical(1))]
    inter <- if (ncol(Xc) >= 2) Xc[, 1] * Xc[, 2] else NULL
    Z <- do.call(cbind, c(parts, list(inter)))
  } else {
    Z <- cbind(poly_cov(), phat_shape())
  }

  Z <- scale(Z)
  Z <- Z[, apply(Z, 2, function(z) all(is.finite(z)) && stats::sd(z) > 1e-8), drop = FALSE]
  if (is.null(dim(Z)) || ncol(Z) < 1)
    stop("The covariate-space basis is degenerate for this fit.")
  qz <- qr(Z)
  if (qz$rank < ncol(Z)) Z <- Z[, qz$pivot[seq_len(qz$rank)], drop = FALSE]

  ZtZ  <- crossprod(Z)
  Zr   <- crossprod(Z, r)
  S    <- as.numeric(t(Zr) %*% solve(ZtZ) %*% Zr)
  ZtU  <- crossprod(Z, sqrt(V) * X)
  ZtOZ <- ZtZ - ZtU %*% solve(crossprod(X, V * X)) %*% t(ZtU)
  lam  <- Re(eigen(solve(ZtZ) %*% ZtOZ, only.values = TRUE)$values)
  lam  <- lam[lam > 1e-9]
  if (!length(lam))
    return(data.frame(Test = "Covariate-space Directed EF", Basis = basis,
                      Test_Statistic = S, df = NA_real_, Method = method,
                      p_value = NA_real_, stringsAsFactors = FALSE))
  data.frame(Test = "Covariate-space Directed EF", Basis = basis,
             Test_Statistic = S, df = sum(lam)^2 / sum(lam^2), Method = method,
             p_value = .def_pvalue(S, lam, method), stringsAsFactors = FALSE)
}


#' Goodness-of-fit evidence features for a fitted model
#'
#' @description
#' Builds the evidence vector used by the learned-ensemble goodness-of-fit test:
#' one-sided z-scores \eqn{\Phi^{-1}(1-p)} from a panel of GOF tests plus the
#' covariate-space directed tests. Larger values mean stronger evidence of misfit.
#'
#' @param object A fitted binary logistic \code{\link[stats]{glm}}.
#' @param tests Character vector of \code{\link{run.all.gof}} test names to use as
#'   panel features (default: a fast partition + DEF-family panel).
#' @return A named numeric vector of evidence features.
#' @seealso \code{\link{deploy.gof}}, \code{\link{cdef.gof}}, \code{\link{run.all.gof}}.
#' @importFrom stats qnorm
#' @export
gof.features <- function(object,
                         tests = c("HL", "HL-equalwidth", "Pigeon-Heyse",
                                   "Tsiatis", "Xie", "EF",
                                   "DEF.poly2", "DEF.poly3", "DEF.stukel")) {
  r  <- suppressWarnings(run.all.gof(object))
  pv <- setNames(r$p_value, r$Test)[tests]
  cp <- suppressWarnings(cdef.gof(object, basis = "poly")$p_value)
  cs <- suppressWarnings(cdef.gof(object, basis = "spline")$p_value)
  p  <- c(pv, CDEF.poly = cp, CDEF.spline = cs)
  z  <- stats::qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6), lower.tail = FALSE)
  z[is.na(z)] <- 0
  z
}


#' Deployable learned-ensemble GOF test via parametric bootstrap
#'
#' @description
#' Turns a pre-trained ensemble \code{meta} into a deployable goodness-of-fit
#' test for \emph{any} fitted model: it scores the model, then calibrates the
#' p-value by a per-dataset parametric bootstrap from the fitted model (so no
#' knowledge of the truth or the data-generating design is required). Validity
#' comes from the bootstrap, independent of how \code{meta} was trained.
#'
#' @param object A fitted binary logistic \code{\link[stats]{glm}}.
#' @param meta A pre-trained scorer: either a function \code{f(features)} returning
#'   a scalar misfit score, or an object with a \code{predict} method consuming a
#'   one-row feature \code{matrix}.
#' @param B Number of parametric-bootstrap resamples (default 99).
#' @param feature_fn Function mapping a fitted glm to its feature vector (default
#'   \code{\link{gof.features}}).
#' @return A one-row \code{data.frame} with the score, \code{B}, and the
#'   bootstrap \code{p_value}.
#' @seealso \code{\link{gof.features}}, \code{\link{cdef.gof}}.
#' @importFrom stats fitted formula glm binomial rbinom predict
#' @export
deploy.gof <- function(object, meta, B = 99, feature_fn = gof.features) {
  stopifnot(inherits(object, "glm"))
  score1 <- function(feat) {
    if (is.function(meta)) return(as.numeric(meta(feat)))
    as.numeric(stats::predict(meta, newx = matrix(feat, 1, dimnames = list(NULL, names(feat))),
                              s = "lambda.min", type = "link"))
  }
  s_obs <- score1(feature_fn(object))
  ph    <- stats::fitted(object)
  dat   <- object$model
  resp  <- names(dat)[1]
  fml   <- stats::formula(object)
  s_star <- vapply(seq_len(B), function(b) {
    dat[[resp]] <- stats::rbinom(length(ph), 1, ph)
    fb <- tryCatch(suppressWarnings(stats::glm(fml, data = dat, family = stats::binomial())),
                   error = function(e) NULL)
    if (is.null(fb) || !isTRUE(fb$converged)) return(NA_real_)
    tryCatch(score1(feature_fn(fb)), error = function(e) NA_real_)
  }, numeric(1))
  s_star <- s_star[is.finite(s_star)]
  p <- (1 + sum(s_star >= s_obs)) / (length(s_star) + 1)
  data.frame(Test = "Deployable learned-ensemble GOF", Score = s_obs, B = length(s_star),
             p_value = p, stringsAsFactors = FALSE)
}
