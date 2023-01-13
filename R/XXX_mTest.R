# ## ----------------------------------------------------------------------------
# ## M-estimator test
# ## ----------------------------------------------------------------------------
#
# #' @title Two sample location test based on M-estimators
# #'
# #' @description \code{m_estimator_test} performs a two-sample permutation or randomization test
# #'              based on M-estimators of location and respective variances
# #'
# #' @template x
# #' @template y
# #' @template alternative
# #' @template delta
# #' @template method_m_test
# #' @template psi
# #' @template k_mest
# #' @template n_rep_m_test
# #' @template na_rm
# #' @template var_test
# #'
# #' @return
# #' A list with class "\code{htest}" containing the following components:
# #' \item{statistic}{the value of the test statistic.}
# #' \item{p.value}{the p-value for the test.}
# #' \item{estimate}{the estimated location difference based on the difference of M-estimates.}
# #' \item{null.value}{the specified hypothesized value of the mean difference.}
# #' \item{alternative}{a character string describing the alternative hypothesis.}
# #' \item{method}{a character string indicating what type of test was performed.}
# #' \item{data.name}{a character string giving the names of the data.}
# #'
# #' @details Details on the tests implemented here can be found in Aboukalam (1992).
# #'
# #' @importFrom Rdpack reprompt
# #'
# #' @references
# #' \insertRef{Abo92robu}{robnptests}
# #'
# #' @examples
# #' x <- rnorm(20); y <- rnorm(20)
# #' m_estimator_test(x, y, psi = "huber", method = "randomization", n.rep = 1000)
# #' m_estimator_test(x, y, psi = "bisquare", method = "randomization", n.rep = 1000)
# #'
# #' @export


# m_estimator_test <- function(x, y, alternative = c("two.sided", "greater", "less"),
#                              delta = ifelse(scale_test, 1, 0),
#                              method = c("randomization", "permutation"),
#                              psi = c("huber", "hampel", "bisquare"), k = .Mpsi.tuning.default(psi),
#                              n.rep = 10000, na.rm = FALSE, scale_test = FALSE) {
#
#   dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
#
#   if (!na.rm & (any(is.na(x)) | any(is.na(y)))) {
#     return(NA)
#   } else if (na.rm & (any(is.na(x)) | any(is.na(y)))) {
#     x <- as.numeric(stats::na.omit(x))
#     y <- as.numeric(stats::na.omit(y))
#   }
#
#   ## If necessary: Transformation to test for difference in scale
#   if (scale_test) {
#     x <- log(x^2)
#     y <- log(y^2)
#     delta <- log(delta^2)
#   }
#
#   alternative <- match.arg(alternative)
#   method <- match.arg(method)
#   psi <- match.arg(psi)
#
#   if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
#     stop ("'delta' must be a single number.")
#   }
#
#   if (length(method == 1) & !(method %in% c("permutation", "randomization"))) {
#     stop (" 'method' must be one of 'permutation' and 'randomization'. ")
#   }
#
#   if (method == "randomization") {
#     randomization <- TRUE
#   } else {
#     randomization <- FALSE
#   }
#
#   stats <- m_test_statistic(x, y - delta, psi = psi, k = k)
#   statistic <- stats$statistic
#
#   estimates <- c(stats$estimates[1], m_est(y, psi = psi, k = k, max.it = 1)$est) #
#
#   distribution <- m_est_perm_distribution(x = x, y = y - delta, randomization = randomization,
#                                          n.rep = n.rep, psi = psi, k1 = k)
#   p.value <- calc_perm_p_value(statistic, distribution, m = length(x), n = length(y),
#                                randomization = randomization, n.rep = n.rep, alternative = alternative)
#
#   if (scale_test) {
#     names(estimates) <- c("M-est. of log(x^2)", "M-est. of log(y^2)")
#     names(delta) <- "ratio of variances"
#     delta <- exp(delta)
#   } else {
#     names(estimates) <- c("M-est. of x", "M-est. of y")
#     names(delta) <- "location shift"
#   }
#   names(statistic) <- ifelse(scale_test, "S", "D")
#
#   if (method == "randomization") {
#     method = paste("Randomization test based on the ", paste0(toupper(substring(psi, 1, 1)), substring(psi, 2, nchar(psi))), "M-estimator")
#   } else method = paste("Exact permutation test based on the", psi, "M-estimator")
#
#   res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
#               estimate = estimates, null.value = delta, alternative = alternative,
#               method = method, data.name = dname)
#
#   class(res) <- "htest"
#
#   return(res)
# }
