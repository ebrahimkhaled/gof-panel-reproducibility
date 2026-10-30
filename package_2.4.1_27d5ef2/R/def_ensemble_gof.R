#' Combine Directed GOF Tests into One Decision (Ensemble)
#'
#' @description
#' Combines the three Directed Ebrahim-Farrington (DEF) basis tests
#' (\code{"poly2"}, \code{"poly3"}, \code{"stukel"}) into a single goodness-of-fit
#' decision, so the user does not have to choose a basis. By default the p-values
#' are combined with the Cauchy Combination Test (CCT), which controls the error
#' rate under the strong dependence between tests computed on the same fitted
#' model. The omnibus EF test can optionally be added to the vote.
#'
#' @details
#' Because the component tests are computed on the same fit, their p-values are
#' strongly dependent. The CCT (\code{combine = "cct"}) has an asymptotic
#' standard-Cauchy null whose tail is robust to this dependence, so it needs no
#' calibration. The \code{"minp"} (Sidak) and \code{"fisher"} rules assume
#' independence and are offered for comparison only; under positive dependence
#' \code{"minp"} is conservative and \code{"fisher"} is anti-conservative, so they
#' should be calibrated by simulation before use (not done here).
#'
#' @param object A fitted binary logistic \code{\link[stats]{glm}}, or a binary
#'   (0/1) vector \code{y} (then supply \code{predicted_probs}).
#' @param predicted_probs Numeric predicted probabilities; required when
#'   \code{object} is a \code{y} vector.
#' @param X Optional design matrix, threaded to \code{\link{def.gof}} for the exact
#'   calibration (only used with the \code{y}/\code{predicted_probs} form).
#' @param components Character vector, a subset of \code{c("poly2","poly3","stukel")}.
#'   Default is all three.
#' @param add_ef Logical; if \code{TRUE}, the omnibus EF p-value (\code{\link{ef.gof}})
#'   is appended to the components. Default \code{FALSE}.
#' @param combine One of \code{"cct"} (default), \code{"minp"}, \code{"fisher"}.
#' @param G Integer number of groups passed to \code{def.gof}/\code{ef.gof} (default 10).
#' @param extra_pvalues Optional named numeric vector of additional p-values to
#'   include (e.g. a Tsiatis test computed elsewhere). Default \code{NULL}.
#'
#' @return A one-row \code{data.frame} with columns \code{Test}, \code{Combiner},
#'   \code{Components}, \code{k}, and \code{p_value}.
#'
#' @references
#' Liu, Y. and Xie, J. (2020). Cauchy combination test. \emph{JASA}, 115(529), 393-402.
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' x <- runif(n, -3, 3)
#' y <- rbinom(n, 1, 1 / (1 + exp(-(0.6 * x))))
#' fit <- glm(y ~ x, family = binomial())
#' def.ensemble.gof(fit)                 # CCT of the three DEF bases
#' def.ensemble.gof(fit, add_ef = TRUE)  # add the omnibus EF
#'
#' @seealso \code{\link{def.gof}}, \code{\link{ef.gof}}.
#' @importFrom stats fitted pchisq
#' @export
def.ensemble.gof <- function(object, predicted_probs = NULL, X = NULL,
                             components = c("poly2", "poly3", "stukel"),
                             add_ef = FALSE,
                             combine = c("cct", "minp", "fisher"),
                             G = 10,
                             extra_pvalues = NULL) {

  combine    <- match.arg(combine)
  components <- match.arg(components, c("poly2", "poly3", "stukel"), several.ok = TRUE)

  # component DEF p-values (Satterthwaite default); flexible input threaded through.
  pv <- vapply(components, function(b) {
    def.gof(object, predicted_probs = predicted_probs, X = X,
            G = G, basis = b, method = "satterthwaite")$p_value
  }, numeric(1))
  names(pv) <- components

  if (isTRUE(add_ef)) {
    if (inherits(object, "glm")) {
      yy <- as.numeric(object$y); pp <- as.numeric(stats::fitted(object))
    } else { yy <- as.numeric(object); pp <- as.numeric(predicted_probs) }
    ef <- ef.gof(yy, pp, G = G)$p_value
    pv <- c(pv, EF = ef)
  }

  if (!is.null(extra_pvalues)) {
    if (!is.numeric(extra_pvalues))
      stop("'extra_pvalues' must be a numeric vector of p-values in [0,1].")
    if (is.null(names(extra_pvalues)))
      names(extra_pvalues) <- paste0("extra", seq_along(extra_pvalues))
    pv <- c(pv, extra_pvalues)
  }

  pcomb <- .combine_pvalues(pv, combine)

  data.frame(
    Test       = "DEF ensemble",
    Combiner   = combine,
    Components = paste(names(pv), collapse = "+"),
    k          = length(pv[is.finite(pv)]),
    p_value    = pcomb,
    stringsAsFactors = FALSE
  )
}

#' EDGES: Cauchy-Combination Ensemble of Directed GOF Tests (Alias)
#'
#' @description
#' Alias for \code{\link{def.ensemble.gof}()}; see the EDGES paper. \code{edges.gof()}
#' is the brand name (EDGES = the Cauchy-combination ensemble of the EDGE directed
#' bases) used in the manuscript. It takes exactly the same arguments as
#' \code{\link{def.ensemble.gof}} and returns exactly the same value; the legacy
#' name \code{def.ensemble.gof()} is retained unchanged for back-compatibility.
#'
#' @param ... Arguments passed on to \code{\link{def.ensemble.gof}} (e.g.
#'   \code{object}, \code{predicted_probs}, \code{X}, \code{components},
#'   \code{add_ef}, \code{combine}, \code{G}, \code{extra_pvalues}).
#'
#' @inherit def.ensemble.gof return references
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' x <- runif(500, -3, 3)
#' y <- rbinom(500, 1, plogis(0.6 * x))
#' fit <- glm(y ~ x, family = binomial())
#' edges.gof(fit)                 # identical to def.ensemble.gof(fit)
#'
#' @seealso \code{\link{def.ensemble.gof}} (legacy name), \code{\link{edge.gof}},
#'   \code{\link{def.gof}}.
#' @export
edges.gof <- function(...) def.ensemble.gof(...)

# Internal: combine a vector of p-values by the chosen rule.
.combine_pvalues <- function(pv, combine) {
  pv <- pv[is.finite(pv)]
  k  <- length(pv)
  if (k == 0) return(NA_real_)
  p  <- pmin(pmax(pv, 1e-15), 1 - 1e-15)
  if (combine == "cct") {
    Tstat <- mean(tan((0.5 - p) * pi))      # Cauchy combination statistic
    0.5 - atan(Tstat) / pi
  } else if (combine == "minp") {
    1 - (1 - min(p))^k                       # Sidak
  } else {                                   # "fisher"
    stats::pchisq(-2 * sum(log(p)), df = 2 * k, lower.tail = FALSE)
  }
}
