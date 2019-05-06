## ----------------------------------------------------------------------------
## Two-sample Hodges-Lehmann test
## ----------------------------------------------------------------------------

#' Two-sample location tests based on two-sample Hodges-Lehmann estimator.
#'
#' @description
#' \code{hl2_test} performs a two-sample location test based on
#' the two-sample Hodges-Lehmann estimator for shift.
#'
#' @param x a numeric vector of observations.
#' @param y a numeric vector of observations.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "\code{two.sided}" (default),
#'                    "\code{greater}" or "\code{less}".
#' @param delta a numeric value indicating the true difference in the location parameter, default is \code{delta = 0}.
#' @param method a character string specifying the test method used, \code{"asymptotic"} for an asymptotic test based on the
#'               normal distribution, \code{"exact"} for an exact and \code{"sampled"} for a randomized permutation test. The exact
#'               permutation test uses all data splits into two samples while the randomized test draws \code{n.rep} random splits
#'               with replacement.
#' @param scale a character string specifying the scale estimator used for standardization of the test statistic,
#'              must be one of \code{"S1"} and \code{"S2"}. Default is \code{"S1"}.
#'              Ignored if \code{method = "asymptotic"}. See details for definition of the scale estimators.
#'
#' @param n.rep an integer value specifying the number of random splits used to calculate
#'              the permutation distribution if \code{method = "sampled"},
#'              ignored if \code{method = "exact"} or \code{method = "asymptotic"}. Default is \code{n.rep = 10000}.
#' @param na.rm a logical value indicating whether NA values in \code{x} and \code{y} should be stripped before the computation proceeds.
#'
#' @details
#' When computing a randomization distribution based on randomly drawn splits with replacement, the results of
#' Smyth & Phipson (2010) to calculate the p-value are used. The test statistics and the asymptotic distribution are taken from Fried & Dehling (2011).
#'
#'  The test statistics for the exact and sampled version of the test is standardized using a robust scale estimator.
#' \code{scale = "S1"} represents use of
#'
#' \deqn{S = med(|X_i - X_j|: 1 \le i < j \le m, |Y_i - Y_j|, 1 <= i < j <= n),}
#'
#' \code{scale = "S2"} uses
#'
#' \deqn{S = med(|Z_i - Z_j|: 1 \le i < j \le m+n) }
#'
#' where \eqn{ Z = ( X_1 - med(X),...,X_m - med(X), Y_1 - med(Y),...,Y_n - med(Y) )'}
#' is the median corrected sample. For more details see Fried & Dehling (2011).
#'
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the estimated difference in means.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of test was performed.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{SmyPhi10perm}{robTests}
#'
#' \insertRef{FriDeh11robu}{robTests}
#'
#' @examples
#' x <- rnorm(20); y <- rnorm(20)
#' hl2_test(x, y, method = "asymptotic", scale = "S1")
#' hl2_test(x, y, method = "sampled", n.rep = 1000, scale = "S2")
#'
#' @export

hl2_test <- function(x, y, alternative = c("two.sided", "greater", "less"),
                     delta = 0, method = c("asymptotic", "exact", "sampled"),
                     scale = c("S1", "S2"), n.rep = 10000,  na.rm = FALSE) {

  stopifnot(is.numeric(x),
            is.numeric(y))

  alternative <- match.arg(alternative)

  scale <- match.arg(scale)

  if (scale == "S1") {
    type <- "D2S1"
  } else if (scale == "S2") {
    type <- "D2S2"
  } else stop(" 'scale' must one of 'S1' and 'S2' ")


  ## Error handling
  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
    stop ("'delta' must be a single number.")
  }

  if (length(method) > 1 & identical(method, c("asymptotic", "exact", "sampled"))) {
    if (length(x) >= 30 & length(y) >= 30) method <- "asymptotic"
    else method <- "sampled"
  }

  if (!(method %in% c("asymptotic", "exact", "sampled"))) {
    stop (" 'method' must be one of 'asymptotic', 'exact' or 'sampled' ")
  }


  if (method %in% c("exact", "sampled")) {
    ## Exact HL2-test using permutation distribution

    ## Results of rob_perm_statistic
    perm.stats <- rob_perm_statistic(x, y - delta, type = type, na.rm = na.rm)

    statistic <- perm.stats$statistic
    # estimates <- perm.stats$estimates

    estimates <- hodges_lehmann_2sample(x, y)

    ## Calculate permutation distribution
    if (method == "sampled") sampled <- TRUE else sampled <- FALSE

    distribution <- perm_distribution(x = x, y = y - delta, type = type, sampled = sampled, n.rep = n.rep)

    ## p-value
    p.value <- calc_perm_p_value(statistic, distribution, m = length(x), n = length(y), sampled = sampled, n.rep = n.rep, alternative = alternative)

    } else if (method == "asymptotic") {

    m <- length(x)
    n <- length(y)
    lambda <- m/(m + n)

    ## Estimation of density at zero for pairwise differences
    xcomb <- utils::combn(x, 2)
    ycomb <- utils::combn(y - delta, 2)

    pwdiffs <- c(xcomb[2, ] - xcomb[1, ], ycomb[2, ] - ycomb[1, ])

    dens <- stats::density(pwdiffs)
    int <- stats::approxfun(dens)(0)

    est <- hodges_lehmann_2sample(x, y - delta)
    statistic <- sqrt(12 * lambda * (1 - lambda)) * int * sqrt(m + n) * est

    estimates <- hodges_lehmann_2sample(x, y)

    p.value <- switch (alternative,
                       two.sided = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
                       greater = stats::pnorm(statistic, lower.tail = FALSE),
                       less = stats::pnorm(statistic, lower.tail = TRUE)
    )
  }

  ## Assign names to results

  names(delta) <- "location shift"
  names(estimates) <- c("HL2 of x and y")
  names(statistic) <- "D"

  if (method == "sampled") {
    method = "Randomization test based on the Two-Sample Hodges-Lehmann estimator"
  } else if (method == "exact") {
    method = "Exact permutation test based on the Two-Sample Hodges-Lehmann estimator"
  } else method = "Asymptotic test based on the Two-Sample Hodges-Lehmann estimator"


  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
              estimate = estimates, null.value = delta, alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}
