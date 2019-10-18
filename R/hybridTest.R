#' @title Hybrid permutation test for location differences
#'
#' @description \code{hybrid_test()} performs a hybrid test based on the p-values of different test statistics
#'
#' @inheritParams hl2_test
#' @param type a character string specifying the type of hybrid test used, must be one of \code{"min1"}, \code{"min2"} and
#' \code{"min3"}. See details for information on the three types.
#' @param k tuning parameter for Huber's M-estimator, default is \code{k = 1.8}, only needed if \code{type \%in\% c("min2", "min3")}
#' @param n.rep an integer value specifying the number of random permutations used to calculate
#'              the permutation distribution of the minimum p-value; default is \code{n.rep = 1000}.
#' @template var_test
#'
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the minimum p-value of the test statistics used.}
#' \item{p.value}{the p-value for the test.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of test was performed.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @details The tests implemented here are introduced in Weichert & Hothorn (2002). They are based on
#' the minimum p-values from different test statistics. The test statistics used are specified in the
#' \code{type}-argument. We have three different types:
#' \describe{
#' \item{min1}{test based on the t-statistic and the 10\% and 20\% trimmed t-statistics}
#' \item{min2}{test based on the t-statistic and a test statistic based on Huber's M-estimator}
#' \item{min3}{test based on the t-statistic, the 20\% trimmed t-statistic and the test based on Huber's M-estimator}
#' }
#'
#' The test statistic is the minumum p-value of the different test statistics.
#' The permutation distribution of the minimum p-value is achieved using the
#' permutation principle according to Efron & Tibshirani (1998).
#'
#' @examples
#' x <- rnorm(10); y <- rnorm(10)
#' hybrid_test(x, y, type = "min1")
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{WeiHot02robu}{robTests}
#'
#' \insertRef{EfrTib98intr}{robTests}
#'
#' @seealso
#'  \code{\link[stats]{t.test}}
#'  \code{\link[robTests]{huber_test}}
#'
#' @export

hybrid_test <- function(x, y, type = c("min1", "min2", "min3"),
                        alternative = c("two.sided", "greater", "less"), delta = 0,
                        k = 1.8, na.rm = FALSE, n.rep = 1000, var.test = FALSE) {

  if (!na.rm & (any(is.na(x)) | any(is.na(y)))) {
    return(NA)
  } else if (na.rm & (any(is.na(x)) | any(is.na(y)))) {
    x <- as.numeric(stats::na.omit(x))
    y <- as.numeric(stats::na.omit(y))
  }

  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
    stop ("'delta' must be a single number.")
  }

  ## If necessary: Transformation to test for difference in scale
  if (var.test) {
    x <- log(x^2)
    y <- log(y^2)
    delta <- log(delta^2)
  }

  type <- match.arg(type)
  alternative <- match.arg(alternative)


  res <- switch(type,
                min1 = min_t_test(x = x, y = y, alternative = alternative, delta = delta,
                                  na.rm = na.rm, n.rep = n.rep),
                min2 = min_c_test(x = x, y = y, alternative = alternative, delta = delta,
                                  k = k, na.rm = na.rm, n.rep = n.rep),
                min3 = min_tc_test(x = x, y = y, alternative = alternative, delta = delta,
                                   k = k, na.rm = na.rm, n.rep = n.rep))

  return(res)
}
