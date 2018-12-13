## ----------------------------------------------------------------------------
## Huber two-sample test
## ----------------------------------------------------------------------------

#' @title Huber 2-sample test
#' @description Performs a 2-sample location test based on Huber's M-estimate
#' @inheritParams hl2_test
#' @param k Tuning parameter for Huber-M-estimate.
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom for the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the Huber M-estimates of \code{x} and \code{y}.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of trimmed t-test was performed.}
#' \item{data.name}{a character string giving the names of the data.}
#' @details 
#' The test is introduced in the context of the hybrid tests in Weichert & Hothorn (2002).
#' 
#' @examples 
#' x <- rnorm(20); y <- rnorm(20)
#' huber_test(x, y)
#' 
#' @seealso 
#' \code{\link[robTests]{min_tc_test}}
#' \code{\link[robTests]{min_t_test}}
#' 
#' @references \insertRef{WeiHot02robu}{robTests}
#'  
#' @export 

huber_test <- function(x, y, delta = 0, k = 1.8, alternative = c("two.sided", "greater", "less")) {
  
  alternative <- match.arg(alternative)
  
  m <- length(x)
  n <- length(y)
  
  hub <- huber_2sample(x, y - delta, k = k)
  
  m.x <- hub$mu.x
  m.y <- hub$mu.y
  
  stand <- hub$s
  
  #statistic <- (m.x$mu - m.y$mu)/sqrt((1/n.x + 1/n.y) * ((n.x - 1) * m.x$s^2 + (n.y - 1) * m.y$s^2)/(n.x + n.y - 2))
  statistic <- (m.x - m.y)/(sqrt(1/m + 1/n) * stand)
  
  if (delta != 0) { 
    estimates <- c(m.x, huber_2sample(x, y, k = k)$mu.y)
  } else estimates <- c(m.x, m.y)
  
  p.value <- switch (alternative,
                     two.sided = 2 * stats::pt(abs(statistic), df = m + n - 2, lower.tail = FALSE),
                     greater = stats::pt(statistic, df = m + n - 2, lower.tail = FALSE),
                     less = stats::pt(statistic, df = m + n - 2, lower.tail = TRUE)
  )
  
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  names(estimates) <- c("Huber-M of x", "Huber-M of y") 
  names(delta) <- "location shift"
  names(statistic) <- "D"
  method <- "Huber-Test"
  parameter <- m + n - 2
  names(parameter) <- "df"
  
  res <- list(statistic = statistic, parameter = parameter, p.value = p.value, 
              estimate = estimates, null.value = delta, alternative = alternative,
              method = method, data.name = dname)
  
  class(res) <- "htest"
  
  return(res)
}



