## ----------------------------------------------------------------------------
## Minimum-t-Test
## ----------------------------------------------------------------------------

#' @title Hybrid permutation test for location difference based on the t- and trimmed t-test
#'
#' @description \code{min_t_test()} performs a hybrid test based on the t-
#'              and Huber-statistic. The significance level is achieved by permutation.
#'
#' @inheritParams hl2_test
#' @param n.rep an integer value specifying the number of random permutations used to calculate
#'              the permutation distribution of the minimum p-value; default is \code{n.rep = 1000}.
#'
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the minimum p-value of t- and Huber-test.}
#' \item{p.value}{the p-value for the test.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of test was performed.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @details The test is introduced in Weichert & Hothorn (2002) and uses the minumum p-value of the t-, the 10%-
#' and the 20%- trimmed t-test as a test statistic.
#' The permutation distribution of the minimum p-value is achieved using the permutation principle according to ???.
#'
#' @examples
#' x <- rnorm(10); y <- rnorm(10)
#' min_t_test(x, y)
#'
#' @references
#' \insertRef{WeiHot02robu}{robTests}
#' @seealso
#'  \code{\link[stats]{t.test}}
#'  \code{\link[robTests]{trimmed_test}}
#' @export
#' @importFrom stats t.test


min_t_test <- function(x, y, alternative = c("two.sided", "greater", "less"), delta = 0,
                       na.rm = FALSE, n.rep = 1000) {

  alternative <- match.arg(alternative)

  z <- c(x, y)
  m <- length(x)
  n <- length(y)

  ## Observed test statistics theta_dach
  t.stat <- stats::t.test(x = x, y = y, mu = delta, alternative = alternative, paired = FALSE, var.equal = TRUE)$statistic
  t10.stat <- trimmed_test(x = x, y = y, delta = delta, gamma = 0.1, alternative = alternative)$statistic
  t20.stat <- trimmed_test(x = x, y = y, delta = delta, gamma = 0.2, alternative = alternative)$statistic

  ## Permutation samples and according test statistics theta_dach_stern
  perm.samples <- replicate(n.rep, sample(z))

  perm.t <- apply(perm.samples, 2, function(z) stats::t.test(x = z[1:m], y = z[(m + 1):(m + n)], mu = delta, alternative = alternative, var.equal = TRUE)$statistic)
  perm.t10 <- apply(perm.samples, 2, function(z) trimmed_test(x = z[1:m], y = z[(m + 1):(m + n)], delta = delta, gamma = 0.1, alternative = alternative)$statistic)
  perm.t20 <- apply(perm.samples, 2, function(z) trimmed_test(x = z[1:m], y = z[(m + 1):(m + n)], delta = delta, gamma = 0.2, alternative = alternative)$statistic)

  pt.obs <- mean(abs(perm.t) >= abs(t.stat))
  pt10.obs <- mean(abs(perm.t10) >= abs(t10.stat))
  pt20.obs <- mean(abs(perm.t20) >= abs(t20.stat))

  obs <- min(pt.obs, pt10.obs, pt20.obs) ### was wollen wir als Statistik ausgeben?
                                         ### und was als estimates?
  statistic <- obs

  pb.t <- sapply(perm.t, function(x) mean(abs(perm.t) >= abs(x)))
  pb.t10 <- sapply(perm.t10, function(x) mean(abs(perm.t10) >= abs(x)))
  pb.t20 <- sapply(perm.t20, function(x) mean(abs(perm.t20) >= abs(x)))

  min.b <- apply(cbind(pb.t, pb.t10, pb.t20), 1, min)

  ## Sachen für die Ausgabe:
  p.value <- mean(min.b <= obs)

  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  names(delta) <- "location shift"
  names(statistic) <- "p"
  method <- "Minimum permutation test based on (trimmed) t-statistics"

  res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
              null.value = delta, alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}
