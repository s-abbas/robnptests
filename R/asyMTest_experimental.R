## ----------------------------------------------------------------------------
## Two-sample location test based on difference of M-estimators
## ----------------------------------------------------------------------------

#' @title Two sample location test based on M-estimators
#'
#' @description \code{m_estimator_test} performs a two-sample location test
#'              based on M-estimators.
#'
#' @template x
#' @template y
#' @template alternative
#' @template delta
#' @template method
#' @template psi
#' @template k_mest
#' @template n_rep
#' @template na_rm
#' @template var_test
#' @template wobble
#' @template seed
#'
#' @details
#' The test statistic is the difference of M-estimates for both samples
#' divided by a pooled tau-estimate for the within-sample variance.
#'
#' When computing a randomization distribution based on randomly drawn splits
#' with replacement, the results of Smyth & Phipson (2010) to calculate the p-value
#' are used. The psi and rho functions, which are needed to obtain the M-estimates,
#' are computed via the implementations in the package \code{\link[robustbase]{robustbase}}.
#' The tau scale estimate is computed with the default parameter settings
#' of the function \code{\link[robustbase]{scaleTau2}}.
#'
#' The test statistic is the difference of the M-estimates for both samples scaled by a pooled estimate for the standard deviation.
#' This estimate is based on the tau-scale estimator. For more details, see the vignette.
#'
#' The distribution of the test statistic is approximated by a standard normal distribution.
#' However, this assumption is only justified under the normality assumption. In case of a non-normal
#' distribution, the test might not keep the desired significance level. The test keeps the
#' level under severe distributions as long as the variance exists. However, under
#' skewed distributions, it tends to be conservative.
#' The test statistic can be corrected by a factor which has to be determined
#' individually for a specific distribution.
#'
#' For \code{var.test = TRUE}, the test compares the two samples for a difference in scale.
#' This is achieved by log-transforming the original observations so that a potential
#' scale difference appears as a location difference between the transformed samples;
#' see Fried (2012). The sample cannot contain zeros due to the necessary log-transformation.
#' If it contains zeros, uniform noise is added to all variables in order to remove zeros.
#' A warning is printed.
#'
#' If the sample has been modified (either because of 0's for \code{var.test = TRUE}, or
#' as \code{wobble = TRUE}, the modified samples can be retrieved using
#'
#' \code{set.seed(wobble.seed); wobble(x, y)}
#'
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom for the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the Huber M-estimates of \code{x} and \code{y}.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of test was performed.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @examples
#' ## Generate random samples
#' set.seed(108)
#' x <- rnorm(20); y <- rnorm(20)
#'
#' ## Asymptotic test based on Huber M-estimator
#' m_test(x, y, method = "asymptotic", psi = "huber")
#'
#' ## Randomization test based on Hampel M-estimator with 1000 random permutations
#' drawn with replacement
#'
#' \dontrun{
#' m_test(x, y, method = "randomization", n.rep = 1000, psi = "hampel")
#' }
#'
#' @export

m_test <- function(x, y, alternative = c("two.sided", "greater", "less"),
                   delta = ifelse(var.test, 1, 0),
                   method = c("asymptotic", "permutation", "randomization"),
                   psi = c("huber", "hampel", "bisquare"),
                   k = robustbase::.Mpsi.tuning.default(psi),
                   n.rep = 10000, na.rm = FALSE,
                   var.test = FALSE, wobble = FALSE, wobble.seed = NULL) {

  ## ___________________________________________________________________________
  ## Error messages
  ## ___________________________________________________________________________

  ## alternative
  if (!all((alternative %in% c("two.sided", "greater", "less")))) {
    stop (" 'alternative' must be one of 'two.sided', 'greater' or 'less'. ")
  }

  ## delta
  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
    stop ("'delta' must be a single number.")
  }

  ## method
  if (!all((method %in% c("asymptotic", "permutation", "randomization")))) {
    stop (" 'method' must be one of 'asymptotic', 'permutation' or 'randomization'. ")
  }

  ## psi
  if (!all((psi %in% c("huber", "hampel", "bisquare")))) {
    stop (" 'psi' must be one of 'huber', 'hampel' or 'bisquare'. ")
  }

  ## ___________________________________________________________________________
  ## Extract default argument values
  ## ___________________________________________________________________________

  alternative <- match.arg(alternative)
  psi <- match.arg(psi)

  ## ___________________________________________________________________________
  ## Check and modify samples
  ## ___________________________________________________________________________

  ## NA handling
  if (!na.rm & (any(is.na(x)) | any(is.na(y)))) {
    return(NA)
  } else if (na.rm & (any(is.na(x)) | any(is.na(y)))) {
    x <- as.numeric(stats::na.omit(x))
    y <- as.numeric(stats::na.omit(y))
  }

  ## Check sample sizes
  if (length(x) < 5 || length(y) < 5) {
    stop("Both samples need at least 5 non-missing values.")
  }

  ## Add random noise if user wants to wobble
  if (wobble) {

    if (is.null(wobble.seed)) wobble.seed <- sample(1e6, 1)
    set.seed(wobble.seed)

    xy <- wobble(x, y)
    x <- xy$x
    y <- xy$y

    warning(paste0("Added random noise to x and y. The seed is ",
                   wobble.seed, "."))
  }

  ## If necessary: Transformation to test for difference in scale
  if (var.test) {
    if (any(c(x, y) == 0)) {

      if (is.null(wobble.seed)) {
        wobble.seed <- sample(1e6, 1)
      }
      set.seed(wobble.seed)

      xy <- wobble(x, y, check = FALSE)
      x <- xy$x
      y <- xy$y

      warning(paste0("Added random noise before log transformation due to zeros in the sample. The seed is ",
                     wobble.seed, "."))
    }
    x <- log(x^2)
    y <- log(y^2)
    delta <- log(delta^2)
  }

  ## ___________________________________________________________________________
  ## Set method for computation of p-value if not specified by user
  ## ___________________________________________________________________________

  ## The method is automatically selected based on the sample sizes
  if (length(method) > 1 & identical(method, c("asymptotic", "permutation", "randomization"))) {
    if (length(x) >= 30 & length(y) >= 30) {
      method <- "asymptotic"
    }
    else {
      method <- "randomization"
      n.rep <- min(choose(length(x) + length(y), length(x)), n.rep)
    }
  }

  ## ___________________________________________________________________________
  ## Test decision
  ## ___________________________________________________________________________

  ## Test statistic and location estimates for both samples
  stats <- m_test_statistic(x, y - delta, psi = psi, k = k)
  statistic <- stats$statistic
  estimates <- stats$estimates

  if (method %in% c("permutation", "randomization")) {
    ## _________________________________________________________________________
    ## Test decision for permutation and randomization test
    ## _________________________________________________________________________

    ## Permutation or randomization distribution
    distribution <- mest_perm_distribution(x = x, y = y - delta, randomization = (method == "randomization"),
                                           n.rep = n.rep, psi = psi, k = k)

    ## p-value
    p.value <- calc_perm_p_value(statistic, distribution, m = length(x), n = length(y),
                                 randomization = (method == "randomization"), n.rep = n.rep, alternative = alternative)
  } else if (method == "asymptotic") {
    ## _________________________________________________________________________
    ## Test decision for asymptotic test
    ## _________________________________________________________________________

    ## p-value
    p.value <- switch (alternative,
                       two.sided = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
                       greater = stats::pnorm(statistic, lower.tail = FALSE),
                       less = stats::pnorm(statistic, lower.tail = TRUE)
    )
  }

  ## ___________________________________________________________________________
  ## Specify output
  ## ___________________________________________________________________________

  ## Names of data sets
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  ## Names of estimates
  if (var.test) {
    names(estimates) <- c("M-est. of log(x^2)", "M-est. of log(y^2)")
    names(delta) <- "ratio of variances"
    delta <- exp(delta)
  } else {
    names(estimates) <- c("M-est. of x", "M-est. of y")
    names(delta) <- "location shift"
  }
  names(statistic) <- ifelse(var.test, "S", "D")

  ## Name of method to compute p-value
  if (method == "randomization") {
    method = paste("Randomization test based on the ", paste0(toupper(substring(psi, 1, 1)), substring(psi, 2, nchar(psi))), "M-estimator")
  } else if (method == "permutation") {
    method = paste("Exact permutation test based on the", psi, "M-estimator")
  } else method = paste("Asymptotic test based on the", psi, "M-estimator")

  ## Output
  res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
              estimate = estimates, null.value = delta, alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}

