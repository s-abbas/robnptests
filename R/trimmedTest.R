## ----------------------------------------------------------------------------
## Trimmed t-Test (Yuen's t-test)
## ----------------------------------------------------------------------------

#' @title Two-sample Trimmed t-test (Yuen's t-Test)
#'
#' @description
#' \code{trimmed_test} performs the two-sample Yuen t-test.
#'
#' @template x
#' @template y
#' @template gamma_trimmed_test
#' @template alternative
#' @template delta
#' @template na_rm
#'
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom for the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the trimmed means of \code{x} and \code{y}.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of test was performed.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @examples
#' ## Generate random samples
#' x <- rnorm(20); y <- rnorm(20)
#'
#' ## Trimmed t-test
#' trimmed_test(x, y, gamma = 0.1)
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{YueDix73appr}{robTests}
#'
#' \insertRef{Yue74trim}{robTests}
#'
#' @export

trimmed_test <- function(x, y, gamma = 0.2,
                         alternative = c("two.sided", "less", "greater"),
                         delta = ifelse(var.test, 1, 0),
                         na.rm = FALSE, var.test = FALSE) {
  ## Error handling
  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
    stop ("'delta' must be a single number.")
  }

  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  ## NA handling
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

  if (var.test) {
    names(estimates) <- c("Trimmed mean of log(x^2)", "Trimmed mean of log(y^2)")
    names(delta) <- "ratio of variances"
    delta <- exp(delta)
  } else {
    names(estimates) <- c("Trimmed mean of x", "Trimmed mean of y")
    names(delta) <- "location shift"
  }

  names(statistic) <- "trimmed t"
  names(df) <- "df"

  method <- "Trimmed two-sample t-test"


  res <- list(statistic = statistic, parameter = df, p.value = p.value,
              estimate = estimates, null.value = delta, alternative = alternative,
              method = method, data.name = dname)
  class(res) <- "htest"

  return(res)
}
