#' Directed Ebrahim-Farrington (DEF) Goodness-of-Fit Test
#'
#' @description
#' Performs the Directed Ebrahim-Farrington (DEF) goodness-of-fit test for a
#' fitted binary logistic regression model. DEF concentrates its power on a small
#' set of calibration-curve "shape" directions by projecting the grouped
#' standardized residuals onto a low-dimensional basis and testing the squared
#' length of that projection.
#'
#' \strong{Naming note:} this test is published under the name \strong{EDGE}
#' (Ebrahim Directed Goodness-of-fit Evaluation), and \code{\link{edge.gof}} is
#' the primary interface going forward. \code{def.gof()} is retained, unchanged,
#' as a fully supported legacy name.
#'
#' @details
#' The observations are sorted by predicted probability and split into \code{G}
#' equal-frequency groups; the standardized grouped residual vector \eqn{r} is
#' projected onto a basis matrix \eqn{Z} of smooth shapes, giving
#' \eqn{S = (Z'r)'(Z'Z)^{-1}(Z'r)}. Its null distribution is a weighted sum of
#' \eqn{\chi^2_1} variables with weights equal to the eigenvalues of
#' \eqn{(Z'Z)^{-1}Z'\Omega Z}, where \eqn{\Omega = I - U(X'WX)^{-1}U'} is the
#' estimation-adjusted covariance of the grouped residuals. The p-value uses a
#' Satterthwaite scaled-\eqn{\chi^2} approximation (default) or Imhof's method
#' (if the \pkg{CompQuadForm} package is installed). Bases: \code{"poly2"},
#' \code{"poly3"} (default), \code{"stukel"}; \code{"ensemble"} runs all three and
#' combines them via \code{\link{def.ensemble.gof}}.
#'
#' @param object A fitted binary logistic \code{\link[stats]{glm}}, or a binary
#'   (0/1) response vector \code{y} (then supply \code{predicted_probs}).
#' @param predicted_probs Numeric predicted probabilities; required when
#'   \code{object} is a \code{y} vector, ignored when it is a glm.
#' @param X Optional design matrix, used only with the \code{y}/\code{predicted_probs}
#'   form: it enables the exact estimation-adjusted (\eqn{\Omega}) calibration
#'   (logit working weights assumed). Without it the conservative \eqn{\chi^2_k}
#'   reference is used and a warning is issued. Ignored when \code{object} is a glm.
#' @param G Integer number of equal-frequency groups (default 10; must be >= 3).
#' @param basis One of \code{"poly3"} (default), \code{"poly2"}, \code{"stukel"},
#'   or \code{"ensemble"}.
#' @param method One of \code{"satterthwaite"} (default) or \code{"imhof"}.
#'
#' @return A one-row \code{data.frame} with columns \code{Test}, \code{Basis},
#'   \code{Test_Statistic} (the statistic \eqn{S}), \code{df}, \code{Method}, and
#'   \code{p_value}. When \code{basis = "ensemble"}, the return is that of
#'   \code{\link{def.ensemble.gof}}.
#'
#' @references
#' Ebrahim, K. E. and El-Kotory, A. Omnibus versus Directed Goodness-of-Fit Tests
#' for Sparse Data in Binary Logistic Regression (companion paper).
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' x <- runif(n, -3, 3)
#' y <- rbinom(n, 1, 1 / (1 + exp(-(0.6 * x))))
#' fit <- glm(y ~ x, family = binomial())
#' def.gof(fit)                       # default poly3 basis
#' def.gof(fit, basis = "stukel")     # tail-shape basis
#' def.gof(fit, basis = "ensemble")   # combine all three (CCT)
#'
#' @seealso \code{\link{ef.gof}}, \code{\link{def.ensemble.gof}}.
#' @importFrom stats fitted predict model.matrix qlogis poly pchisq
#' @export
def.gof <- function(object, predicted_probs = NULL, X = NULL, G = 10,
                    basis  = c("poly3", "poly2", "stukel", "ensemble"),
                    method = c("satterthwaite", "imhof")) {

  basis  <- match.arg(basis)
  method <- match.arg(method)
  if (!is.numeric(G) || length(G) != 1 || G < 3) {
    stop("'G' must be a single integer >= 3.")
  }

  if (basis == "ensemble")
    return(def.ensemble.gof(object, predicted_probs = predicted_probs, X = X, G = G))

  # --- accept either a fitted glm, OR (y, predicted_probs[, X]) ---
  if (inherits(object, "glm")) {
    if (object$family$family != "binomial")
      stop("'object' must be a binomial glm (or pass y, predicted_probs, X).")
    y   <- as.numeric(object$y)
    ph  <- pmin(pmax(as.numeric(stats::fitted(object)), 1e-6), 1 - 1e-6)
    eta <- as.numeric(stats::predict(object, type = "link"))
    dmu <- object$family$mu.eta(eta)
    X   <- stats::model.matrix(object)
    naive <- FALSE
  } else {
    if (!is.numeric(object))
      stop("'object' must be a fitted binomial glm or a numeric (0/1) y vector.")
    y <- as.numeric(object)
    if (is.null(predicted_probs))
      stop("Provide 'predicted_probs' when 'object' is not a glm.")
    ph  <- pmin(pmax(as.numeric(predicted_probs), 1e-6), 1 - 1e-6)
    dmu <- ph * (1 - ph)
    naive <- is.null(X)
    if (naive)
      warning("def.gof: no model/design matrix supplied; using the conservative ",
              "chi-square reference (Omega = I). Pass the fitted glm or X for the ",
              "exact estimation-adjusted test.")
  }
  n <- length(y)
  if (!all(y %in% c(0, 1))) stop("DEF needs a binary (0/1) response.")
  if (length(ph) != n) stop("'object' (y) and 'predicted_probs' lengths differ.")
  if (G > n) stop("'G' cannot exceed the number of observations.")

  V <- ph * (1 - ph)
  w <- dmu^2 / V

  # --- equal-frequency groups by fitted probability ---
  grp  <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  idx  <- split(seq_len(n), grp)
  Gn   <- length(idx)
  og   <- vapply(idx, function(I) sum(y[I]),   numeric(1))
  eg   <- vapply(idx, function(I) sum(ph[I]),  numeric(1))
  Vg   <- vapply(idx, function(I) sum(V[I]),   numeric(1))
  pbar <- vapply(idx, function(I) mean(ph[I]), numeric(1))
  r    <- (og - eg) / sqrt(Vg)

  # --- estimation-adjusted covariance ---
  if (naive) {
    Omega <- diag(Gn)
  } else {
    U     <- t(vapply(idx, function(I) colSums(dmu[I] * X[I, , drop = FALSE]),
                      numeric(ncol(X)))) / sqrt(Vg)
    Omega <- diag(Gn) - U %*% solve(crossprod(X, w * X)) %*% t(U)
  }

  # --- shape basis Z ---
  Z <- .def_basis(pbar, basis)
  Z <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]
  if (ncol(Z) < 1)
    stop("The chosen basis is degenerate for this fit. Try basis = 'poly3' or a larger G.")

  ZtZ <- crossprod(Z)
  Zr  <- crossprod(Z, r)
  S   <- as.numeric(t(Zr) %*% solve(ZtZ) %*% Zr)

  lam <- Re(eigen(solve(ZtZ) %*% (t(Z) %*% Omega %*% Z), only.values = TRUE)$values)
  lam <- lam[lam > 1e-9]
  if (length(lam) == 0) {
    return(data.frame(Test = "Directed Ebrahim-Farrington", Basis = basis,
                      Test_Statistic = S, df = NA_real_, Method = method,
                      p_value = NA_real_, stringsAsFactors = FALSE))
  }

  pval <- .def_pvalue(S, lam, method)
  nu   <- sum(lam)^2 / sum(lam^2)

  data.frame(Test = "Directed Ebrahim-Farrington", Basis = basis,
             Test_Statistic = S, df = nu, Method = method, p_value = pval,
             stringsAsFactors = FALSE)
}

# Internal: build the shape-basis matrix Z from the group mean probabilities.
.def_basis <- function(pbar, basis) {
  if (basis %in% c("poly2", "poly3")) {
    deg <- if (basis == "poly2") 2L else 3L
    if (length(unique(round(pbar, 8))) < deg + 1)
      stop("Too few distinct group probabilities for basis = '", basis,
           "'. Use a smaller polynomial basis or a larger sample/G.")
    as.matrix(stats::poly(pbar, deg))
  } else {
    e <- stats::qlogis(pbar)
    cbind(e, e^2 * (e >= 0), -e^2 * (e < 0))
  }
}

# Internal: p-value of S under sum_j lambda_j chi^2_1.
.def_pvalue <- function(S, lam, method) {
  if (method == "imhof") {
    if (!requireNamespace("CompQuadForm", quietly = TRUE))
      stop("method = 'imhof' requires the 'CompQuadForm' package. ",
           "Install it, or use method = 'satterthwaite'.")
    p <- CompQuadForm::imhof(S, lam)$Qq
    return(min(max(p, 0), 1))
  }
  cc <- sum(lam^2) / sum(lam)
  nu <- sum(lam)^2 / sum(lam^2)
  stats::pchisq(S / cc, df = nu, lower.tail = FALSE)
}
