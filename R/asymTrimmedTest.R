## ----------------------------------------------------------------------------
## Asymmetrically trimmed tests
## ----------------------------------------------------------------------------

#' Two-sample location tests based on asymmetrically trimmed means
#'
#' @description
#' \code{asym_trimmed_test()} performs a two-sample location test by using
#' asymmetrically trimmed means based on different skewness-
#' selector statistics for both samples.
#'
#'
#' @template x
#' @template y
#' @template type
#' @template alternative
#' @template delta
#' @template method
#' @template n_rep
#' @template na_rm
#'
#'
#' @details
#' The test statistic is an analogue to the one of the ordinary t-test, where
#' the difference of the sample means is replaced by the difference of asymmetrically trimmed means
#' and the pooled empirical standard deviation is replaced by a pooled winsorized standard deviation
#' which uses the asymmetrically trimmed means.
#'
#' The number of observations trimmed from the lower and the upper end of the sample
#' depends on the skewness. Reed & Stark (2004) suggest three possible skewness-selector statistics.
#' The argument \code{type} specifies which one is used.
#'
#' \describe{
#' \item{"Q2"}{This skewness-selector statistic compares the difference between the upper and the lower 5\% of
#'             the observations to the difference between the upper and lower 50\%.}
#' \item{"SK2"}{This skewness-selector statistic compares the difference between the sample minimum and
#'              median to the difference between the median and the maximum.}
#' \item{"SK5"}{This skewness-selector statistic compares the difference between the sample minimum and
#'              mean to the difference between the mean and the maximum.}
#' }
#'
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
#' \insertRef{ReeSta04robu}{robTests}
#'
#' @examples
#' set.seed(108)
#' x <- rnorm(20); y <- rnorm(20)
#' asym_trimmed_test(x, y, type = "Q2", method = "asymptotic")
#' asym_trimmed_test(x, y, type = "SK2", method = "asymptotic")
#' asym_trimmed_test(x, y, type = "SK5", method = "asymptotic")
#'
#' asym_trimmed_test(x, y, type = "SK5", method = "sampled")
#'
#' @export

asym_trimmed_test <- function(x, y, type = c("Q2", "SK2", "SK5"),
                              alternative = c("two.sided", "greater", "less"),
                              delta = 0,
                              method = c("asymptotic", "exact", "sampled"),
                              n.rep = 10000, na.rm = FALSE) {

  alternative <- match.arg(alternative)
  method <- match.arg(method)
  type <- match.arg(type)

  ## Error handling
  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
    stop ("'delta' must be a single number.")
  }

  # if (length(method) > 1 & identical(method, c("asymptotic", "exact", "sampled"))) {
  #   if (length(x) >= 30 & length(y) >= 30) method <- "asymptotic"
  #   else method <- "sampled"
  # }

  if (!(method %in% c("asymptotic", "exact", "sampled"))) {
    stop (" 'method' must be one of 'asymptotic', 'exact' or 'sampled' ")
  }

  ## Test statistic
  t.stat <- asym_trimmed_t(x, y - delta, type = type, na.rm = na.rm)

  statistic <- t.stat$statistic
  estimates <- t.stat$estimates
  df <- t.stat$df

  if (method %in% c("exact", "sampled")) {
    ## Calculate permutation distribution
    if (method == "sampled") sampled <- TRUE else sampled <- FALSE

    distribution <- asym_trimmed_perm_distribution(x = x, y = y - delta, type = type,
                                      sampled = sampled, n.rep = n.rep)

    ## p-value
    p.value <- calc_perm_p_value(statistic, distribution, m = length(x), n = length(y),
                                 sampled = sampled, n.rep = n.rep, alternative = alternative)

  } else if (method == "asymptotic") {
      p.value <- switch (alternative,
                         two.sided = 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE),
                         greater = stats::pt(statistic, df = df, lower.tail = FALSE),
                         less = stats::pt(statistic, df = df, lower.tail = TRUE)
    )
  }

  ## Assign names to results
  names(estimates) <- c("Trimmed mean of x", "Trimmed mean of y")
  names(delta) <- "location shift"
  names(statistic) <- "D"

  if (method == "sampled") {
    method = paste("Randomization test based on the", type, "selector statistic")
  } else if (method == "exact") {
    method = paste("Exact permutation test based on the", type, "selector statistic")
  } else method = paste("Asymptotic test based on the", type, "selector statistic")


  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
              estimate = estimates, null.value = delta, alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}
