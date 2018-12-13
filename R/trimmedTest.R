## ----------------------------------------------------------------------------
## Trimmed t-Test (Yuen's t-test)
## ----------------------------------------------------------------------------

#' @title Two-sample Trimmed t-test (Yuen's t-Test)
#' 
#' @description 
#' \code{trimmed_test} performs the two-sample Yuen t-test.
#' 
#' @param x numeric vector of observations.
#' @param y numeric vector of observations.
#' @param gamma numeric value in [0, 0.5] specifying the fraction of observations to be trimmed from each end of the sample before calculating the mean. Values of trim outside that range are taken as the nearest endpoint.
#' @param alternative character string specifying the alternative hypothesis, must be one of "\code{two.sided}" (default), "\code{greater}" or "\code{less}".
#' @param delta numeric indicating the true difference in means
#' @param na.rm a logical value indicating whether NA values in \code{x} should be stripped before the computation proceeds.
#' 
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom for the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the estimated difference in means.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{data.name}{a character string giving the names of the data.}
#' 
#' @examples 
#' x <- rnorm(20); y <- rnorm(20)
#' trimmed_test(x, y, gamma = 0.1)
#' 
#' @references
#' \insertRef{YueDix73appr}{robTests}
#' 
#' \insertRef{Yue74trim}{robTests}
#' 
#' @export
#' 
# trimmed_test <- function(x, ...) {
#   UseMethod("trimmed_test")
# }
# 
# #' @rdname trimmed_test
# #' @method trimmed_test default
# #' 
# #' @export

trimmed_test <- function(x, y, gamma = 0.2, alternative = c("two.sided", "less", "greater"), 
                                 delta = 0, na.rm = FALSE) {
  ## Error handling
  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
    stop ("'delta' must be a single number.")
  }
  names(delta) <- "difference in means" ### BB: Ist das hier überflüssig?
  
  ## NA handling 
  if (!na.rm & any(is.na(x))) {
    return(NA)
  } else if (na.rm & any(is.na(x))) {
    x <- stats::na.omit(x)
  }
  
  alternative <- match.arg(alternative)
  
  ## Results of trimmed_t
  t.stats <- trimmed_t(x, y, delta = delta, gamma = gamma, na.rm = na.rm)
  
  statistic <- t.stats$statistic
  estimates <- t.stats$estimates
  df <- t.stats$df
  
  ## p-value
  p.value <- switch (alternative,
                     two.sided = 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE),
                     greater = stats::pt(statistic, df = df, lower.tail = FALSE),
                     less = stats::pt(statistic, df = df, lower.tail = TRUE)
  )
  
  ## Assign names to results
  names(statistic) <- "trimmed t"
  names(estimates) <- c("Trimmed mean of x", "Trimmed mean of y")
  names(delta) <- "difference in means"
  names(df) <- "df"

  method <- "Trimmed two-sample t-test"

  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  res <- list(statistic = statistic, parameter = df, p.value = p.value, 
              estimate = estimates, null.value = delta, alternative = alternative, 
              method = method, data.name = dname)
  class(res) <- "htest"
  
  return(res)
}
