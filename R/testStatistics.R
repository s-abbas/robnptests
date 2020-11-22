## ----------------------------------------------------------------------------
## Calculation of the different test statistics
## ----------------------------------------------------------------------------

#' @title Test statistic for the two-sample Yuen t-test
#'
#' @description
#' \code{trimmed_t} calculates the test statistic of the two-sample Yuen test.
#'
#' @template x
#' @template y
#' @template gamma_trimmed_test
#' @template delta
#' @template na_rm
#'
#' @return
#' A list containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{estimates}{the trimmed means for both samples.}
#' \item{df}{the degrees of freedom for the test statistic.}
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{YueDix73appr}{robTests}
#'
#' \insertRef{Yue74trim}{robTests}
#'
#' @export

trimmed_t <- function(x, y, gamma = 0.2, delta = 0, na.rm = FALSE) {

  ## Trimmed means
  x.trim <- trim_mean(x, gamma = gamma, na.rm = na.rm)
  y.trim <- trim_mean(y - delta, gamma = gamma, na.rm = na.rm)
  estimates <- c(x.trim, y.trim)

  ## Scale estimator
  var.x <- win_var(x, gamma = gamma, na.rm = na.rm)
  var.y <- win_var(y, gamma = gamma, na.rm = na.rm)

  h.x <- var.x$h
  h.y <- var.y$h
  var.x <- var.x$var
  var.y <- var.y$var

  m <- length(x)
  n <- length(y)
  pool.var <- ((m - 1) * var.x + (n - 1) * var.y)/(h.x + h.y - 2)

  ## Degrees of freedom
  df <- h.x + h.y - 2

  ## Test statistic
  statistic <- (x.trim - y.trim) / sqrt(pool.var * (1/h.x + 1/h.y))

  res <- list(statistic = statistic, estimates = estimates, df = df)

  return(res)
}



#' @title Robust permutation statistics based on robust location estimators
#'
#' @description \code{rob_perm_statistic()} calculates test statistics for robust
#' permutation/randomization tests based on the sample median, the one-sample
#' Hodges-Lehmann estimator, or the two-sample Hodges-Lehmann estimator.
#'
#' @template x
#' @template y
#' @template type_rob_perm
#' @template na_rm
#'
#' @return A list containing the following components:
#'         \item{statistic}{the selected test statistic.}
#'         \item{abs.statistic}{the absolute value of the selected test statistic.}
#'         \item{estimates}{estimate of location for each sample if available.}
#'
#' @details The test statistics returned by \code{rob_perm_statistic} are of the form
#'          \deqn{D_i/S_j} where the D_i, i = 1,...,3, are different estimators of location and the S_j, j = 1,...,4 are
#'          estimates for the mutual sample scale. See \insertCite{FriDeh11robu;textual}{robTests}
#'          or the vignette (\code{vignette(robTests-vignette)}) for details.
#'
#' @examples
#' ## Generate random samples
#' set.seed(108)
#' x <- rnorm(20); y <- rnorm(20)
#'
#' ## Compute HL21 statistic
#' #rob_perm_statistic(x, y, type = "HL21")
#'
#' @references
#' \insertRef{FriDeh11robu}{robTests}
#'
#' @export

rob_perm_statistic <- function(x, y,
                          type = c("HL11", "HL12", "HL21", "HL22", "MED1", "MED2"),
                          na.rm = FALSE) {
  type <- match.arg(type)

  if(!(type %in% c("HL11", "HL12", "HL21", "HL22", "MED1", "MED2"))) {
    stop("type needs to be one of 'HL11', 'HL12', 'HL21', 'HL22', 'MED1', 'MED2'.")
  }

  switch(type,
         HL11 = {
            est.x <- hodges_lehmann(x, na.rm = na.rm)
            est.y <- hodges_lehmann(y, na.rm = na.rm)
            loc <- est.x - est.y
            sd <- rob_var(x, y, na.rm = na.rm, type = "S1")
            res <- loc/sd
         },
         HL12 = {
            est.x <- hodges_lehmann(x, na.rm = na.rm)
            est.y <- hodges_lehmann(y, na.rm = na.rm)
            loc <- est.x - est.y
            sd <- rob_var(x, y, na.rm = na.rm, type = "S2")
            res <- loc/sd
         },
         HL21 = {
            est.x <- est.y <- NULL
            loc <- hodges_lehmann_2sample(x, y, na.rm = na.rm)
            sd <- rob_var(x, y, na.rm = na.rm, type = "S1")
            res <- loc/sd
         },
         HL22 = {
            est.x <- est.y <- NULL
            loc <- hodges_lehmann_2sample(x, y, na.rm = na.rm)
            sd <- rob_var(x, y, na.rm = na.rm, type = "S2")
            res <- loc/sd
         },
         MED1 = {
            est.x <- stats::median(x, na.rm = na.rm)
            est.y <- stats::median(y, na.rm = na.rm)
            loc <- est.x - est.y
            sd <- rob_var(x, y, na.rm = na.rm, type = "S3")
            res <- loc/sd
        },
         MED2 = {
            est.x <- stats::median(x, na.rm = na.rm)
            est.y <- stats::median(y, na.rm = na.rm)
            loc <- est.x - est.y
            sd <- rob_var(x, y, na.rm = na.rm, type = "S4")
            res <- loc/sd
      })

  return(list(statistic = res, abs.statistic = abs(res), estimates = c(est.x, est.y)))
}


#' @title Calculate the M-test-statistic
#'
#' @description \code{m_test_statistic} calculates the test statistics for tests based on M-estimators.
#'
#' @template x
#' @template y
#' @template psi
#' @template k_mest
#' @template scaleTau2
#'
#' @return A list containing the following components:
#'         \item{statistic}{standardized test statistic.}
#'         \item{estimates}{M-estimates of location for both \code{x} and \code{y}.}
#'
#' @export


m_test_statistic <- function(x, y, psi, k = robustbase::.Mpsi.tuning.default(psi),
                             ...) {
  ## Sample sizes
  m <- length(x)
  n <- length(y)

  ## M-estimates for both samples
  est.x <- m_est(x = x, psi = psi, k = k)$est
  est.y <- m_est(x = y, psi = psi, k = k)$est

  ## Estimator for \nu
  psi.x <- robustbase::Mpsi((x - est.x)/stats::mad(x), psi = psi, cc = k)
  rho.x <- robustbase::Mpsi((x - est.x)/stats::mad(x), psi = psi, cc = k, deriv = 1)

  psi.y <- robustbase::Mpsi((y - est.y)/stats::mad(y), psi = psi, cc = k)
  rho.y <- robustbase::Mpsi((y - est.y)/stats::mad(y), psi = psi, cc = k, deriv = 1)

  nu.x <- mean(psi.x^2)/(mean(rho.x)^2)
  nu.y <- mean(psi.y^2)/(mean(rho.y)^2)

  ## Test statistic
  return(list(statistic = (est.x - est.y) / sqrt((n * robustbase::scaleTau2(x, consistency = TRUE, ...)^2 * nu.x + m * robustbase::scaleTau2(y, consistency = TRUE, ...)^2 * nu.y) / (m * n)),
              estimates = c(est.x, est.y)))
}


