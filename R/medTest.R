## ----------------------------------------------------------------------------
## Two sample test based on median differences
## ----------------------------------------------------------------------------

#' Two-sample tests based on the difference of the one-sample Hodges-Lehmann
#' estimator.
#' 
#' @description 
#' \code{med_test} performs a two-sample test for the location shift based on
#' the difference of the sample medians for both samples.
#' 
#' @param scale a character string specifying the scale estimator used for standardization of the test statistic,
#'              must be one of \code{"S3"} and \code{"S4"}. Default is \code{"S3"}. 
#'              Ignored if \code{method = "asymptotic"}. See details for definition of the scale estimators.
#'             
#' @inheritParams hl2_test
#' 
#' @details
#' If the number of replications to determine the permutation distribution is drawn with replacement, we use the results of
#' Smyth & Phipson (2010) to calculate the p-value.
#' The test statistics and the asymptotic distribution are taken from Fried & Dehling (2011).
#' 
#' The test statistics for the exact and sampled version of the test is standardized using a robust scale estimator. 
#' \code{scale = "S3"} represents use of
#' 
#' \deqn{S = 2 * ( |X_1 - med(X)|,...,|X_m - med(X)|, |Y_1 - med(Y)|,...,|Y_n - med(Y)| ),} 
#'
#' \code{scale = "S4"} uses
#' 
#' \deqn{S = 1.4826 * ( med( |X_1 - med(X)|,...,|X_m - med(X)| ) + med( |Y_1 - med(Y)|,...,|Y_n - med(Y)| ). }
#' 
#' For more details see Fried & Dehling (2011). 
#' 
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the sample medians of \code{x} and \code{y}.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of trimmed t-test was performed.}
#' \item{data.name}{a character string giving the names of the data.}
#' 
#' @references
#' \insertRef{SmyPhi10perm}{robTests}
#' 
#' \insertRef{FriDeh11robu}{robTests}
#' 
#' @examples
#' x <- rnorm(20); y <- rnorm(20)
#' med_test(x, y, method = "asymptotic", scale = "S3")
#' med_test(x, y, method = "sampled", n.rep = 1000, scale = "S4")
#' 
#' @export

med_test <- function(x, y, alternative = c("two.sided", "greater", "less"), delta = 0,
                     method = c("asymptotic", "exact", "sampled"), 
                     scale = c("S3", "S4"), n.rep = 10000, 
                     na.rm = FALSE) {
 
  alternative <- match.arg(alternative)
  scale <- match.arg(scale)
  
  if (scale == "S3") {
    type <- "D3S3"
  } else if (scale == "S4") {
    type <- "D3S4"
    } else stop(" 'scale' must one of 'S3' and 'S4' ")
  
  ## Error handling
  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
    stop ("'delta' must be a single number.")
  }
  
  if (length(method) > 1 & identical(method, c("asymptotic", "exact", "sampled"))) {
    if (length(x) >= 30 & length(y) >= 30) method <- "asymptotic"
    else method <- "sampled"
  }
  
  if (!(method %in% c("asymptotic", "exact", "sampled"))) {
    stop (" 'method' must be one of 'asymptotic', 'exact' or 'sampled' ")
  }
  
  if (method %in% c("exact", "sampled")) {
    
    ## Results of rob_perm_statistic
    perm.stats <- rob_perm_statistic(x, y - delta, type = type, na.rm = na.rm)
    
    statistic <- perm.stats$statistic
    estimates <- perm.stats$estimates
    
    if (delta != 0) estimates[2] <- stats::median(y)
    
    ## Calculate permutation distribution
    
    if (method == "sampled") sampled <- TRUE else sampled <- FALSE
    
    distribution <- perm_distribution(x = x, y = y - delta, type = type, sampled = sampled, 
                                      n.rep = n.rep)
    
    ## p-value
    p.value <- calc_perm_p_value(statistic, distribution, m = length(x), 
                                 n = length(y), sampled = sampled, n.rep = n.rep,
                                 alternative = alternative)
  
  } else if (method == "asymptotic") { 
    
    med.x <- stats::median(x, na.rm = na.rm)
    med.y <- stats::median(y - delta, na.rm = na.rm)
    
    diff <- c(x - med.x, y - delta - med.y) 
    
    dens <- stats::approxfun(stats::density(diff))
    med <- dens(0)  
    
    m <- length(x)
    n <- length(y)
    
    if (delta != 0) estimates <- c(med.x, stats::median(y)) else estimates <- c(med.x, med.y)
    est <- med.y - med.x
    
    statistic <- sqrt(m*n/(m+n)) * 2 * med * est
    
    p.value <- switch (alternative,
                       two.sided = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
                       greater = stats::pnorm(statistic, lower.tail = FALSE),
                       less = stats::pnorm(statistic, lower.tail = TRUE)
    )
   
  }
  
  ## Assign names to results
  
  names(estimates) <- c("Median of x", "Median of y")
  names(delta) <- "location shift"
  names(statistic) <- "D"
  
  if (method == "sampled") {
    method = "Randomization test based on sample medians"
  } else if (method == "exact") {
    method = "Exact permutation test based on sample medians"
  } else method = "Asymptotic test based on sample medians"
  
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  res <- list(statistic = statistic, parameter = NULL, p.value = p.value, 
              estimate = estimates, null.value = delta, alternative = alternative, 
              method = method, data.name = dname)
  
  class(res) <- "htest"
  
  return(res)
}
