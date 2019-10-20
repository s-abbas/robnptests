## ----------------------------------------------------------------------------
## One-sample Hodges-Lehmann test
## ----------------------------------------------------------------------------

#' Two-sample location tests based on one-sample Hodges-Lehmann estimator
#'
#' @description
#' \code{hl1_test} performs a two-sample location test based on
#' the difference of the one-sample Hodges-Lehmann estimators for both samples.
#'
#' @template x
#' @template y
#' @template alternative
#' @template delta
#' @template method
#' @template scale_hl
#' @template n_rep
#' @template na_rm
#' @template var_test
#'
#' @details
#' When computing a randomization distribution based on randomly drawn splits with replacement, the results of
#' Smyth & Phipson (2010) to calculate the p-value are used. The test statistics and the asymptotic distribution are taken from Fried & Dehling (2011).
#'
#' The test statistics for the exact and sampled version of the test is standardized using a robust scale estimator.
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
#' \item{estimate}{the estimated sample means of \code{x} and \code{y}.}
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
#' @import utils
#'
#' @examples
#' x <- rnorm(20); y <- rnorm(20)
#' hl1_test(x, y, method = "asymptotic", scale = "S1")
#' hl1_test(x, y, method = "sampled", n.rep = 1000, scale = "S2")
#'
#' @export

hl1_test <- function(x, y, alternative = c("two.sided", "greater", "less"), delta = 0,
                     method = c("asymptotic", "exact", "sampled"), scale = c("S1", "S2"),
                     n.rep = 10000, na.rm = FALSE,
                     var.test = FALSE) {

  alternative <- match.arg(alternative)

  scale <- match.arg(scale)

  if (!na.rm & (any(is.na(x)) | any(is.na(y)))) {
    return(NA)
  } else if (na.rm & (any(is.na(x)) | any(is.na(y)))) {
    x <- as.numeric(stats::na.omit(x))
    y <- as.numeric(stats::na.omit(y))
  }

  ## If necessary: Transformation to test for difference in scale
  if (var.test) {
    x <- log(x^2)
    y <- log(y^2)
    delta <- log(delta^2)
  }

  if (scale == "S1") {
    type <- "D1S1"
  } else if (scale == "S2") {
    type <- "D1S2"
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
    ## Results of rob_perm_statistic
    perm.stats <- rob_perm_statistic(x, y - delta, type = type, na.rm = na.rm)

    statistic <- perm.stats$statistic
    estimates <- perm.stats$estimates
    if (delta != 0) estimates[2] <- hodges_lehmann(y)

    ## Calculate permutation distribution
    if (method == "sampled") sampled <- TRUE else sampled <- FALSE

    distribution <- perm_distribution(x = x, y = y - delta, type = type,
                                      sampled = sampled, n.rep = n.rep)

    ## p-value
    p.value <- calc_perm_p_value(statistic, distribution, m = length(x), n = length(y),
                                 sampled = sampled, n.rep = n.rep, alternative = alternative)

  } else if (method == "asymptotic") {

    m <- length(x)
    n <- length(y)

    # pairwise differences the density estimate is calculated from:
    xcomb <- utils::combn(x, 2)
    ycomb <- utils::combn(y - delta, 2)
    pwdiffs <- c(xcomb[2, ] - xcomb[1, ], ycomb[2, ] - ycomb[1, ])
    dens <- stats::density(pwdiffs)
    dens <- stats::approxfun(dens)

    int <- dens(0)

    estimates <- c(hodges_lehmann(x), hodges_lehmann(y - delta))
    statistic <- sqrt(12*m*n/(m+n)) * int * (estimates[2] - estimates[1])

    if (delta != 0) estimates[2] <- hodges_lehmann(y)

    p.value <- switch (alternative,
                       two.sided = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
                       greater = stats::pnorm(statistic, lower.tail = FALSE),
                       less = stats::pnorm(statistic, lower.tail = TRUE)
    )
  }
  ## Assign names to results

  names(estimates) <- c("HL1 of x", "HL1 of y")
  names(delta) <- "location shift"
  names(statistic) <- "D"

  if (method == "sampled") {
    method = "Randomization test based on the Hodges-Lehmann estimator"
  } else if (method == "exact") {
    method = "Exact permutation test based on the Hodges-Lehmann estimator"
  } else method = "Asymptotic test based on the Hodges-Lehmann estimator"


  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
              estimate = estimates, null.value = delta, alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}
