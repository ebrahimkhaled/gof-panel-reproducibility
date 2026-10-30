#' EDGE: Directed Goodness-of-Fit Test for Binary Logistic Regression
#'
#' @description
#' \code{edge.gof()} is the primary interface to the EDGE test (Ebrahim Directed
#' Goodness-of-fit Evaluation): a grouped, \emph{directed} goodness-of-fit test
#' for binary logistic regression under sparse data. EDGE projects the grouped
#' standardized residuals onto a small pre-specified basis of calibration shapes
#' (cubic \code{"poly3"} by default) and refers the resulting quadratic form to
#' its closed-form weighted chi-squared null distribution -- no refit, no
#' resampling, no tuning.
#'
#' \code{edge.gof()} computes exactly the same statistic as the legacy name
#' \code{\link{def.gof}} (retained for backward compatibility); the returned
#' \code{Test} label is \code{"EDGE"}.
#'
#' @inheritParams def.gof
#'
#' @return A one-row \code{data.frame} with columns \code{Test} (\code{"EDGE"}),
#'   \code{Basis}, \code{Test_Statistic}, \code{df}, \code{Method}, and
#'   \code{p_value}, as documented in \code{\link{def.gof}}.
#'
#' @references
#' Ebrahim, E. K. and El-Kotory, A. (2026). EDGE: a closed-form directed
#' goodness-of-fit test for sparse logistic regression. Manuscript.
#'
#' Ebrahim, E. K. and El-Kotory, A. (2026). A directional Hosmer-Lemeshow
#' goodness-of-fit test for sparse logistic regression. arXiv:2607.15454.
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' x <- runif(500, -3, 3)
#' y <- rbinom(500, 1, plogis(0.6 * x))
#' fit <- glm(y ~ x, family = binomial())
#' edge.gof(fit)                      # default cubic basis, G = 10
#' edge.gof(fit, basis = "stukel")    # Stukel-shape basis
#'
#' @seealso \code{\link{def.gof}} (legacy name), \code{\link{ef.gof}},
#'   \code{\link{def.ensemble.gof}}, \code{\link{run.all.gof}}.
#' @export
edge.gof <- function(object, predicted_probs = NULL, X = NULL, G = 10,
                     basis = "poly3", method = "satterthwaite") {
  out <- def.gof(object, predicted_probs = predicted_probs, X = X, G = G,
                 basis = basis, method = method)
  out$Test <- "EDGE"
  out
}
