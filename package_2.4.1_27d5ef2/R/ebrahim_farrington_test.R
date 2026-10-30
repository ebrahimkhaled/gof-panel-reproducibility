#' Ebrahim-Farrington Goodness-of-Fit Test for Logistic Regression
#'
#' @description
#' Performs the Ebrahim-Farrington goodness-of-fit test for logistic regression models.
#' This test is particularly effective for binary data and sparse datasets, providing
#' an improved alternative to the traditional Hosmer-Lemeshow test.
#'
#' @details
#' The Ebrahim-Farrington test is based on Farrington's (1996) theoretical framework
#' but simplified for practical implementation with binary data. The test uses a
#' modified Pearson chi-square statistic with data-dependent grouping, where
#' observations are grouped by their predicted probabilities.
#'
#' For binary data (when \code{G} is specified), the test automatically groups
#' observations into \code{G} groups based on predicted probabilities and applies
#' the simplified Ebrahim-Farrington statistic:
#'
#' \deqn{Z_{EF} = \frac{T_{EF} - (G - 2)}{\sqrt{2(G-2)}}}
#'
#' where \eqn{T_{EF}} is the modified Pearson chi-square statistic, and \eqn{G}
#' is the number of groups.
#'
#' For grouped data (when \code{m} is provided), the test applies the original
#' Farrington test with full variance calculations.
#'
#' @param y A fitted binary logistic \code{glm} (then \code{predicted_probs} is
#'   taken from it automatically), or a numeric vector of binary responses (0/1)
#'   for binary data / counts of successes for grouped data.
#' @param predicted_probs Numeric vector of predicted probabilities from the
#'   logistic regression model. Must be same length as \code{y}.
#' @param model Optional \code{glm} object. Required only for the original
#'   Farrington test with grouped data (when \code{m} is provided and \code{G} is NULL).
#' @param m Optional numeric vector of trial counts for each observation
#'   (for grouped data). If NULL, data is assumed to be binary.
#' @param G Optional integer specifying the number of groups for binary data
#'   grouping. Default is 10. If NULL, no grouping is performed and \code{m}
#'   must be provided.
#' @param method Reference distribution for the grouped EF statistic:
#'   \code{"chisq"} (default) refers \eqn{T_{EF}} to a \eqn{\chi^2_{G-2}}
#'   distribution; \code{"normal"} uses the standardized \eqn{Z_{EF}} (the
#'   behaviour of package versions <= 1.0.0).
#'
#' @return
#' A data frame with the following columns:
#' \item{Test}{Character string identifying the test performed}
#' \item{Test_Statistic}{Numeric value of the standardized test statistic}
#' \item{p_value}{Numeric p-value for the test}
#'
#' @note
#' \itemize{
#'   \item For binary data with automatic grouping (\code{G} specified): Use the
#'         Ebrahim-Farrington test which is computationally efficient and doesn't
#'         require the model specification.
#'   \item For grouped data (\code{m} provided): Use the original Farrington test
#'         which requires the fitted model object.
#'   \item The test statistic follows a standard normal distribution under the
#'         null hypothesis of adequate model fit.
#'   \item For binary data with \code{m=1} for all observations and no grouping,
#'         the test is not applicable and will return a p-value of 1.
#' }
#'
#' @references
#' Farrington, C. P. (1996). On Assessing Goodness of Fit of Generalized Linear Models to Sparse Data. *Journal of the Royal Statistical Society. Series B (Methodological)*, 58(2), 349-360.
#' Ebrahim, K. E. (2025). Goodness-of-Fits Tests and Calibration Machine Learning Algorithms for Logistic Regression Model with Sparse Data. *Master's Thesis*, Alexandria University.
#' Hosmer, D. W., & Lemeshow, S. (1980). A goodness-of-fit test for the multiple logistic regression model. *Communications in Statistics - Theory and Methods*, 9(10), 1043–1069. https://doi.org/10.1080/03610928008827941
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' # Example 1: Binary data with automatic grouping (Ebrahim-Farrington test)
#' set.seed(123)
#' n <- 500
#' x <- rnorm(n)
#' linpred <- 0.5 + 1.2 * x
#' prob <- 1 / (1 + exp(-linpred))
#' y <- rbinom(n, 1, prob)
#' 
#' # Fit logistic regression
#' model <- glm(y ~ x, family = binomial())
#' predicted_probs <- fitted(model)
#' 
#' # Perform Ebrahim-Farrington test with 10 groups
#' result <- ef.gof(y, predicted_probs, G = 10)
#' print(result)
#' 
#' # Example 2: Compare with different number of groups
#' result_4 <- ef.gof(y, predicted_probs, G = 4)
#' result_20 <- ef.gof(y, predicted_probs, G = 20)
#' 
#' # Example 3: Grouped data (original Farrington test)
#' # Note: This requires actual grouped data with trials > 1
#' \dontrun{
#' # Simulated grouped data
#' n_groups <- 50
#' m_trials <- sample(5:20, n_groups, replace = TRUE)
#' x_grouped <- rnorm(n_groups)
#' linpred_grouped <- -0.5 + 1.0 * x_grouped
#' prob_grouped <- 1 / (1 + exp(-linpred_grouped))
#' y_grouped <- rbinom(n_groups, m_trials, prob_grouped)
#' 
#' # Fit model for grouped data
#' data_grouped <- data.frame(successes = y_grouped, trials = m_trials, x = x_grouped)
#' model_grouped <- glm(cbind(successes, trials - successes) ~ x, 
#'                      data = data_grouped, family = binomial())
#' predicted_probs_grouped <- fitted(model_grouped)
#' 
#' # Original Farrington test
#' result_grouped <- ef.gof(y_grouped, predicted_probs_grouped, 
#'                          model = model_grouped, m = m_trials)
#' print(result_grouped)
#' }
#'
#' @seealso
#' \code{\link[ResourceSelection]{hoslem.test}} for the Hosmer-Lemeshow test
#'
#' @export
ef.gof <- function(y, predicted_probs = NULL, model = NULL, m = NULL, G = 10,
                   method = c("chisq", "normal")) {

  method <- match.arg(method)

  # The first argument may be a fitted binomial glm; extract y and predictions.
  if (inherits(y, "glm")) {
    if (y$family$family != "binomial") stop("ef.gof: the model must be a binomial glm.")
    predicted_probs <- as.numeric(stats::fitted(y))
    y <- as.numeric(y$y)
  }
  if (is.null(predicted_probs)) {
    stop("Provide 'predicted_probs' (or pass a fitted glm as the first argument).")
  }

  # Input validation
  n <- length(y)
  if (n != length(predicted_probs)) {
    stop("Length mismatch: 'y' and 'predicted_probs' must have the same length")
  }
  
  if (!is.numeric(y) || !is.numeric(predicted_probs)) {
    stop("'y' and 'predicted_probs' must be numeric vectors")
  }
  
  if (any(predicted_probs < 0) || any(predicted_probs > 1)) {
    stop("'predicted_probs' must be between 0 and 1")
  }
  
  if (!all(y %in% c(0, 1)) && is.null(m)) {
    stop("For binary data, 'y' must contain only 0s and 1s")
  }
  
  # Determine which test to perform
  if (!is.null(G)) {
    # Ebrahim-Farrington test with automatic grouping
    if (!is.null(model) || !is.null(m)) {
      message("Note: When using Ebrahim-Farrington grouping (G specified), ",
              "the 'model' and 'm' parameters are not needed. ",
              "Ignoring these parameters for grouping mode.")
    }
    
    if (!is.numeric(G) || length(G) != 1 || G < 2) {
      stop("'G' must be a single integer >= 2")
    }
    
    if (G > n) {
      stop("Number of groups 'G' cannot exceed sample size 'n'")
    }
    
    return(.ebrahim_farrington_grouped(y, predicted_probs, G, method))
    
  } else {
    # Original Farrington test for grouped data
    if (is.null(m)) {
      stop("For the original Farrington test without grouping, ",
           "you must provide the 'm' vector (trial counts). ",
           "For binary data, use G parameter for automatic grouping.")
    }
    
    if (is.null(model)) {
      stop("For the original Farrington test with grouped data, ",
           "you must provide the fitted 'model' object.")
    }
    
    return(.farrington_original(y, predicted_probs, model, m))
  }
}

# Internal function for Ebrahim-Farrington test with grouping
.ebrahim_farrington_grouped <- function(y, predicted_probs, G, method = "chisq") {
  n <- length(y)
  
  # Sort data by predicted probabilities
  sorted_indices <- order(predicted_probs)
  sorted_y <- y[sorted_indices]
  sorted_probs <- predicted_probs[sorted_indices]
  
  # Create groups of approximately equal size
  group_size <- floor(n / G)
  remainder <- n %% G
  
  # Create group assignments
  group_assignments <- rep(1:G, each = group_size)
  if (remainder > 0) {
    group_assignments <- c(group_assignments, rep(G, remainder))
  }
  
  # Calculate grouped statistics
  y_grouped <- tapply(sorted_y, group_assignments, sum)
  m_grouped <- tapply(sorted_y, group_assignments, length)
  probs_grouped <- tapply(sorted_probs, group_assignments, mean)
  
  # Convert to numeric vectors
  y_group <- as.numeric(y_grouped)
  m_group <- as.numeric(m_grouped)
  prob_group <- as.numeric(probs_grouped)
  
  # Apply numerical stability bounds
  prob_bounded <- pmax(pmin(prob_group, 0.999999), 0.000001)
  
  # Calculate Ebrahim-Farrington statistic
  variance_terms <- m_group * prob_bounded * (1 - prob_bounded)
  pearson_components <- ((y_group - m_group * prob_bounded)^2) / variance_terms
  farrington_correction <- ((1 - 2 * prob_bounded) / variance_terms) * 
                          (y_group - m_group * prob_bounded)
  
  test_statistic_raw <- sum(pearson_components) - sum(farrington_correction)
  
  # Simplified expectation and variance for Ebrahim-Farrington test
  expected_value <- G - 2
  variance_value <- 2 * (G - 2)
  
  # Handle edge case where variance is too small
  if (variance_value <= 1e-8 || is.na(variance_value)) {
    z_statistic <- 0
    p_value <- 1
    warning("Variance too small for reliable test. Consider using more groups.")
  } else {
    z_statistic <- (test_statistic_raw - expected_value) / sqrt(variance_value)
    # chisq (default): refer T_EF to chi-square_{G-2}; normal: the standardized Z_EF
    p_value <- if (identical(method, "normal")) {
      1 - stats::pnorm(z_statistic)
    } else {
      stats::pchisq(test_statistic_raw, df = G - 2, lower.tail = FALSE)
    }
  }

  # Return results
  data.frame(
    Test = "Ebrahim-Farrington",
    Test_Statistic = z_statistic,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}

# Internal function for original Farrington test
.farrington_original <- function(y, predicted_probs, model, m) {
  n <- length(y)
  
  # Validate inputs
  if (length(m) != n) {
    stop("Length mismatch: 'y' and 'm' must have the same length")
  }
  
  if (!inherits(model, "glm")) {
    stop("'model' must be a glm object")
  }
  
  # Apply numerical stability bounds
  prob_bounded <- pmax(pmin(predicted_probs, 0.999999), 0.000001)
  
  # Get design matrix
  design <- stats::model.matrix(model)
  p <- ncol(design)
  
  # Calculate variance terms
  variance_terms <- m * prob_bounded * (1 - prob_bounded)
  
  # Calculate test statistic components
  pearson_components <- ((y - m * prob_bounded)^2) / variance_terms
  farrington_correction <- ((1 - 2 * prob_bounded) / variance_terms) * 
                          (y - m * prob_bounded)
  
  test_statistic_raw <- sum(pearson_components) - sum(farrington_correction)
  
  # Calculate expected value and variance for original Farrington test
  # This requires the information matrix
  tryCatch({
    w_matrix <- diag(as.vector(variance_terms))
    xpxi <- solve(t(design) %*% w_matrix %*% design)
    varlinpr <- diag(design %*% xpxi %*% t(design))
    
    expected_value <- (n - p) + sum(prob_bounded * (1 - prob_bounded) * varlinpr)
    variance_value <- 2 * (1 - (p / n)) * sum((m - 1) / m)
    
  }, error = function(e) {
    warning("Could not compute full Farrington variance. Using simplified version.")
    expected_value <- n - p
    variance_value <- 2 * (n - p) * mean((m - 1) / m)
  })
  
  # Handle cases where test is not applicable
  if (variance_value <= 1e-8 || is.na(variance_value)) {
    z_statistic <- 0
    p_value <- 1
    if (all(m == 1)) {
      warning("Farrington test is not applicable to purely binary data (all m=1). ",
              "Consider using G parameter for automatic grouping.")
    }
  } else {
    z_statistic <- (test_statistic_raw - expected_value) / sqrt(variance_value)
    p_value <- 1 - stats::pnorm(z_statistic)
  }
  
  # Return results
  data.frame(
    Test = "Farrington-Original",
    Test_Statistic = z_statistic,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
} 