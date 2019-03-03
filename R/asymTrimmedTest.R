## ----------------------------------------------------------------------------
## Asymmetrically trimmed test based on Q2 skewness-selector statistic
## ----------------------------------------------------------------------------

#' Two-sample location tests based on the difference of the Q2 skewness-selector
#' statistic.
#'
#' @description
#' \code{asym_trimmed_test()} performs a two-sample test for the location shift based on
#' the difference of the asymmetrically trimmed means based on different skewness-
#' selector statistic for both samples.
#'
#'
#' @inheritParams hl2_test
#' @param  type specifies the skewness selector statistic used for trimming, must be in
#'        \code{"Q2", "SK2"} and \code{"SK5"}. Default is \code{"Q2"}.
#'
#'
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the estimated sample means of \code{x} and \code{y}.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of trimmed t-test was performed.}
#' \item{data.name}{a character string giving the names of the data.}

#' @references
#' \insertRef{SmyPhi10perm}{robTests}
#' \insertRef{ReeSta04robu}{robTests}

#'
#' @import utils
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
