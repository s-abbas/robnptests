## ----------------------------------------------------------------------------
## Two-sample Hodges-Lehmann test
## ----------------------------------------------------------------------------

#' Two-sample location tests based on two-sample Hodges-Lehmann estimator.
#'
#' @description
#' \code{hl2_test} performs a two-sample location test based on
#' the two-sample Hodges-Lehmann estimator for shift.
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
#'  The test statistics for the permutation and randomization version of the test is standardized using a robust scale estimator.
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
#' \item{estimate}{the estimated location difference based on the two-sample Hodges-Lehmann estimator.}
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
#' hl2_test(x, y, method = "randomization", n.rep = 1000, scale = "S2")
#'
#' @export

hl2_test <- function(x, y, alternative = c("two.sided", "greater", "less"),
                     delta = ifelse(var.test, 1, 0), method = c("asymptotic", "permutation", "randomization"),
                     scale = c("S1", "S2"), n.rep = 10000,  na.rm = FALSE,
                     var.test = FALSE) {

  stopifnot(is.numeric(x),
            is.numeric(y))

  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

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
  method <- match.arg(method)
  scale <- match.arg(scale)



  if (scale == "S1") {
    type <- "HL21"
  } else if (scale == "S2") {
    type <- "HL22"
  } else stop(" 'scale' must one of 'S1' and 'S2' ")


  ## Error handling
  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
    stop ("'delta' must be a single number.")
  }

  if (!all(method %in% c("asymptotic", "permutation", "randomization"))) {
    stop (" 'method' must be one of 'asymptotic', 'permutation' or 'randomization'. ")
  }

  ## If no choice is made regarding the computation of the p-value, the method
  ## is automatically selected based on the sample sizes
  if (length(method) > 1 & identical(method, c("asymptotic", "permutation", "randomization"))) {
    if (length(x) >= 30 & length(y) >= 30) {
      method <- "asymptotic"
    }
    else {
      method <- "randomization"
      n.rep <- min(choose(length(x) + length(y), length(x)), n.rep)
    }
  }

  if (method %in% c("permutation", "randomization")) {
    ## Exact HL2-test using permutation distribution

    ## Results of rob_perm_statistic
    perm.stats <- rob_perm_statistic(x, y - delta, type = type, na.rm = na.rm)

    statistic <- perm.stats$statistic
    # estimates <- perm.stats$estimates

    estimates <- hodges_lehmann_2sample(x, y)

    ## Calculate permutation distribution
    if (method == "randomization") {
      randomization <- TRUE
    } else {
      randomization <- FALSE
    }

    distribution <- perm_distribution(x = x, y = y - delta, type = type, randomization = (method == randomization), n.rep = n.rep)

    ## p-value
    p.value <- calc_perm_p_value(statistic, distribution, m = length(x), n = length(y), randomization = randomization, n.rep = n.rep, alternative = alternative)

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
  if (var.test) {
    names(estimates) <- c("HL2 of log(x^2) and log(y^2)")
    names(delta) <- "ratio of variances"
    delta <- exp(delta)
  } else {
    names(estimates) <- c("HL2 of x and y")
    names(delta) <- "location shift"
  }

  names(statistic) <- ifelse(var.test, "S", "D")

  if (method == "randomization") {
    method = "Randomization test based on the Two-Sample Hodges-Lehmann estimator"
  } else if (method == "permutation") {
    method = "Exact permutation test based on the Two-Sample Hodges-Lehmann estimator"
  } else method = "Asymptotic test based on the Two-Sample Hodges-Lehmann estimator"


  res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
              estimate = estimates, null.value = delta, alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}
