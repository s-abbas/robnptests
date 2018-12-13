## ----------------------------------------------------------------------------
## M-estimator test
## ----------------------------------------------------------------------------

#' @title Two sample location test based on M-estimators
#' 
#' @description \code{m_estimator_test} performs a two-sample permutation or randomization test
#'              based on M-estimators of location and respective variances
#'              
#' @inheritParams hl2_test
#' @param method a character string specifying the test method used: \code{"exact"} for an exact and \code{"sampled"}
#'               for a randomized permutation test. The exact 
#'               permutation test uses all data permutations while the randomized test draws \code{n.rep} random permutations
#'               with replacement. 
#' @param psi kernel used for optimization, must be one of \code{"bisquare"}, \code{"hampel"} and \code{"huber"},
#'            defaults to \code{"huber"}.
#' @param k tuning parameter(s) for the respective kernel function, defaults to parameters implemented in .Mpsi.tuning.default(psi) from 
#'          \code{robustbase} package.
#'          
#' @param n.rep an integer value specifying the number of random permutations used to calculate
#'              the permutation distribution if \code{method = "sampled"}, 
#'              ignored if \code{method = "exact"}. Default is 10000.
#'              
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the estimated difference in means.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of trimmed t-test was performed.}
#' \item{data.name}{a character string giving the names of the data.}
#' 
#' @details Details on the tests implemented here can be found in Aboukalam (1992).
#' @references 
#' \insertRef{Abo92robu}{robTests}
#' 
#' @import robustbase
#' @examples 
#' x <- rnorm(20); y <- rnorm(20)
#' m_estimator_test(x, y, psi = "huber", method = "sampled", n.rep = 1000)
#' m_estimator_test(x, y, psi = "bisquare", method = "sampled", n.rep = 1000)
#' @export 


m_estimator_test <- function(x, y, alternative = c("two.sided", "greater", "less"), delta = 0, 
                             method = c("sampled", "exact"), 
                             psi = c("huber", "hampel", "bisquare"), k = .Mpsi.tuning.default(psi),
                             n.rep = 10000, na.rm = FALSE) {
  
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  psi <- match.arg(psi)
  
  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
    stop ("'delta' must be a single number.")
  }
  
  if (length(method == 1) & !(method %in% c("exact", "sampled"))) {
    stop (" 'method' must be one of 'exact' and 'sampled' ")
  }
  
  if (method == "sampled") sampled <- TRUE else sampled <- FALSE
    
    stats <- m_test_statistic(x, y - delta, psi = psi, k = k)
    statistic <- stats$statistic
    
    estimates <- c(stats$estimates[1], m_est(y, psi = psi, k = k, max.it = 1)$est) # 
    
    distribution <- mest_perm_distribution(x = x, y = y - delta, sampled = sampled, 
                                           n.rep = n.rep, psi = psi, k1 = k)
    p.value <- calc_perm_p_value(statistic, distribution, m = length(x), n = length(y),
                                 sampled = sampled, n.rep = n.rep, alternative = alternative)
    
    names(estimates) <- c("M-est. of x", "M-est. of y") ## y - delta? 
    names(delta) <- "location shift"
    names(statistic) <- "D"
    
    if (method == "sampled") {
      method = paste("Randomization test based on the ", paste0(toupper(substring(psi, 1, 1)), substring(psi, 2, nchar(psi))), "M-estimator") # ist das okay so mit den Kleinbuchstaben?
    } else method = paste("Exact permutation test based on the", psi, "M-estimator")
    
    
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    
    res <- list(statistic = statistic, parameter = NULL, p.value = p.value, 
                estimate = estimates, null.value = delta, alternative = alternative,
                method = method, data.name = dname)
    
    class(res) <- "htest"
  
    return(res)
}
