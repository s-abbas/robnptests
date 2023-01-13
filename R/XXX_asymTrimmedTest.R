#' #' @title Asymmetrically trimmed mean
#' #'
#' #' @description
#' #' \code{asym_trimmed_mean} calculates an asymmetrically trimmed mean of a sample.
#' #' In contrast to an ordinary trimmed mean, the numbers of observations removed
#' #' from the lower and the upper end of the sample do not need to be equal.
#' #'
#' #' @template x
#' #' @template type_skewness
#' #' @template na_rm
#' #'
#' #' @details
#' #' The number of observations trimmed from the lower and the upper end of the sample
#' #' depends on the skewness.
#' #'
#' #' Three skewness-selector statistics suggested by Reed and Stark (1996) have been
#' #' implemented. The argument \code{type} specifies which one is used.
#' #'
#' #' \describe{
#' #' \item{"Q2"}{This skewness-selector statistic compares the difference between the upper and the lower 5\% of
#' #'             the observations to the difference between the upper and lower 50\%.}
#' #' \item{"SK2"}{This skewness-selector statistic compares the difference between the sample minimum and
#' #'              median to the difference between the median and the maximum.}
#' #' \item{"SK5"}{This skewness-selector statistic compares the difference between the sample minimum and
#' #'              mean to the difference between the mean and the maximum.}
#' #' }
#' #'
#' #' @return
#' #' The asymmetrically trimmed mean.
#' #'
#' #' @importFrom Rdpack reprompt
#' #'
#' #' @references
#' #' \insertRef{ReeSta96hing}{robnptests}
#' #'
#' #' @examples
#' #' ## Generate random sample
#' #' set.seed(108)
#' #' x <- rnorm(10)
#' #'
#' #' ## Compute asymmetrically trimmed means
#' #' asym_trimmed_mean(x, type = "Q2")
#' #' asym_trimmed_mean(x, type = "SK2")
#' #' asym_trimmed_mean(x, type = "SK5")
#' #'
#' #' @export
#'
#' asym_trimmed_mean <- function(x, type = c("Q2", "SK2", "SK5"), na.rm = FALSE) {
#'   ## NA handling
#'   if (!na.rm & any(is.na(x))) {
#'     return(NA)
#'   } else if (na.rm & any(is.na(x))) {
#'     x <- as.vector(stats::na.omit(x))
#'   }
#'
#'   type <- match.arg(type)
#'
#'   ## Sample size and sorted sample
#'   m <- length(x)
#'   x.sort <- sort(x)
#'
#'   ## Determine lower trimming proportion
#'   if (type == "Q2") {
#'
#'     ## Do not trim, if sample size is too small
#'     if (floor(0.05 * m) == 0) {
#'       return(mean(x))
#'     }
#'
#'     gamma <- 0.1
#'
#'     U.005 <- mean(x.sort[(m - floor(0.05 * m) + 1):m])
#'     L.005 <- mean(x.sort[1:floor(0.05 * m)])
#'     U.05 <- mean(x.sort[(m - floor(0.5 * m) + 1):m])
#'     L.05 <- mean(x.sort[1:floor(0.5 * m)])
#'
#'     gamma.lower <- gamma * (U.005 - L.005)/(U.005 - L.005 + U.05 - L.05)
#'   } else if (type == "SK2") {
#'     gamma <- 0.1
#'
#'     gamma.lower <- gamma * (min(x) - stats::median(x))/(min(x) - max(x))
#'   } else if (type == "SK5") {
#'     gamma <- 0.25
#'
#'     gamma.lower <- gamma * (min(x) - mean(x))/(min(x) - max(x))
#'   }
#'
#'   ## Determine upper trimming proportion
#'   gamma.upper <- gamma - gamma.lower
#'
#'   ## Calculate number of values to be trimmed from both ends
#'   r.lower <- floor(m * gamma.lower) + 1
#'   r.upper <- m - floor(m * gamma.upper)
#'
#'   ## Asymmetrically trimmed mean
#'   res <- mean(x.sort[(r.lower + 1):r.upper])
#'
#'   return(res)
#' }
#'
#' #' @title Asymmetrically trimmed variance
#' #'
#' #'
#' #' @description
#' #' \code{asym_win_var} calculates the asymmetrically trimmed variance using different
#' #' skewness selector statistics.
#' #'
#' #'
#' #' @template x
#' #' @template type_skewness
#' #' @template na_rm
#' #'
#' #'
#' #' @return A list containing the following items:
#' #' \item{var}{The asymmetrically trimmed variance, and}
#' #' \item{h}{the degrees of freedom used for tests based on  asymmetrically trimmed means and the
#' #' asymmetrically winsorized variance.}
#' #'
#' #'
#' #' @examples
#' #' ## Generate random sample
#' #' set.seed(108)
#' #' x <- rnorm(10)
#' #'
#' #' ## Compute asymmetrical winsorized variance
#' #' asym_win_var(x, type = "SK5")
#' #'
#' #'
#' #' @export
#'
#' asym_win_var <- function(x, type = c("Q2", "SK2", "SK5"), na.rm = FALSE) {
#'   ## NA handling
#'   if (!na.rm & any(is.na(x))) {
#'     return(NA)
#'   } else if (na.rm & any(is.na(x))) {
#'     x <- as.vector(stats::na.omit(x))
#'   }
#'
#'   ## Sample size and ordered sample
#'   m <- length(x)
#'   x.sort <- sort(x)
#'
#'   ## Calculate asymmetrically trimmed mean
#'   mean.x <- asym_trimmed_mean(x, type = type)
#'
#'   if (type == "Q2") {
#'     if (floor(0.05 * m) == 0) {
#'       ## If sample size is too small for trimming, calculate ordinary
#'       ## empirical variance
#'       res <- stats::var(x)
#'
#'       return(list(var = res, h = m))
#'     }
#'
#'     gamma <- 0.1
#'
#'     ## Determine lower trimming proportion
#'     U.005 <- mean(x.sort[(m - floor(0.05 * m) + 1):m])
#'     L.005 <- mean(x.sort[1:floor(0.05 * m)])
#'     U.05 <- mean(x.sort[(m - floor(0.5 * m) + 1):m])
#'     L.05 <- mean(x.sort[1:floor(0.5 * m)])
#'
#'     gamma.lower <- gamma * (U.005 - L.005)/(U.005 - L.005 + U.05 - L.05)
#'   } else if (type == "SK2") {
#'     gamma <- 0.1
#'
#'     gamma.lower <- gamma * (min(x) - stats::median(x))/(min(x) - max(x))
#'   } else if (type == "SK5") {
#'     gamma <- 0.25
#'
#'     gamma.lower <- gamma * (min(x) - mean(x))/(min(x) - max(x))
#'   }
#'
#'   ## Determine upper trimming proportion
#'   gamma.upper <- gamma - gamma.lower
#'
#'   ## Replace observations in sorted sample to compute winsorized variance
#'   m.lower <- 1 + floor(gamma.lower * m + 0.5)
#'   m.upper <- m - floor(gamma.upper * m + 0.5)
#'
#'   x.sort[1:(m.lower - 1)] <- x.sort[m.lower]
#'   x.sort[(m.upper + 1):m] <- x.sort[m.upper]
#'
#'   ## Calculate variance
#'   res <- 1/(m.upper - m.lower + 1) * sum((x.sort - mean.x)^2)
#'
#'   h <- m.upper - m.lower + 1
#'
#'   return(list(var = res, h = h))
#' }
#'
#'
#' #' @title Test statistic for the asymmetrically trimmed test
#' #'
#' #' @description
#' #' \code{asym_trimmed_t} calculates the test statistic of the asymmetrically trimmed test with a given skewness selector statistic
#' #'
#' #' @template x
#' #' @template y
#' #' @template type_skewness
#' #' @template delta
#' #' @template na_rm
#' #'
#' #' @return
#' #' A list containing the following components:
#' #' \item{statistic}{the value of the test statistic.}
#' #' \item{estimates}{the asymmetrically trimmed means for both samples.}
#' #' \item{df}{the degrees of freedom for the test statistic.}
#' #'
#' #' @references
#' #' \insertRef{YueDix73appr}{robnptests}
#' #'
#' #' \insertRef{Yue74trim}{robnptests}
#' #'
#' #' @export
#'
#' asym_trimmed_t <- function(x, y, type, delta = 0, na.rm = FALSE) {
#'
#'   ## Trimmed means
#'   x.trim <- asym_trimmed_mean(x, type = type, na.rm = na.rm)
#'   y.trim <- asym_trimmed_mean(y - delta, type = type, na.rm = na.rm)
#'   estimates <- c(x.trim, y.trim)
#'
#'   ## Scale estimator
#'   scale.x <- asym_win_var(x, type = type, na.rm = na.rm)
#'   scale.y <- asym_win_var(y, type = type, na.rm = na.rm)
#'
#'   h.x <- scale.x$h
#'   h.y <- scale.y$h
#'   scale.x <- scale.x$var
#'   scale.y <- scale.y$var
#'
#'   m <- length(x)
#'   n <- length(y)
#'   pool.var <- ((h.x - 1) * scale.x + (h.y - 1) * scale.y)/(h.x + h.y - 2)
#'
#'   ## Degrees of freedom
#'   df <- h.x + h.y - 2
#'
#'   ## Test statistic
#'   statistic <- (x.trim - y.trim) / sqrt(pool.var * (1/h.x + 1/h.y))
#'
#'   res <- list(statistic = statistic, estimates = estimates, df = df)
#'
#'   return(res)
#' }
#'
#'
#' #' @title Permutation distribution for asymmetrically trimmed statistics
#' #'
#' #' @description \code{asym_trimmed_perm_distribution} calculates the permutation distribution for the asymmetrically
#' #' trimmed statistics from \code{asym_trimmed_test}.
#' #'
#' #' @template x
#' #' @template y
#' #' @template type_skewness
#' #' @template randomization
#' #' @template n_rep
#' #'
#' #' @return Vector with permutation distribution.
#' #'
#' #' @export
#'
#' asym_trimmed_perm_distribution <- function(x, y, type, randomization = FALSE, n.rep = NULL) {
#'   ## Sample sizes
#'   m <- length(x)
#'   n <- length(y)
#'
#'   ## Splits in two samples
#'   if (!randomization) {
#'     complete <- c(x, y)
#'     splits <- gtools::combinations((m + n), m, 1:(m + n))
#'
#'     distribution <- apply(splits, 1, function(s) asym_trimmed_t(x = complete[s], y = complete[-s], type)$statistic)
#'
#'   } else if (randomization) {
#'     splits <- replicate(n.rep, sample(c(x, y)))
#'
#'     distribution <- apply(splits, 2, function(s) asym_trimmed_t(x = s[1:m], y = s[(m + 1):(m + n)], type)$statistic)
#'   }
#'
#'   return(distribution)
#' }
#'
#' #' Calculation of permutation p-value
#' #'
#' #' @description
#' #' \code{calc_perm_p_value} calculates the permutation p-value following Phipson & Smyth (2010).
#' #'
#' #' @template statistic
#' #' @template distribution
#' #' @template m
#' #' @template n
#' #' @template randomization
#' #' @template n_rep
#' #' @template alternative
#' #'
#' #' @return
#' #' p.value for the specified alternative.
#' #'
#' #' @references
#' #' \insertRef{PhiSmy10perm}{robnptests}
#' #'
#' #' @export
#'
#' calc_perm_p_value <- function(statistic, distribution, m, n, randomization, n.rep, alternative) {
#'
#'   ## Number of permutations leading to test statistic at least as extreme
#'   ## as observed
#'   A <- switch(alternative,
#'               two.sided = sum(abs(distribution) >= abs(statistic)),
#'               greater = sum(distribution >= statistic),
#'               less = sum(distribution <= statistic)
#'   )
#'
#'   ## Computation of p-value
#'   if (randomization) {
#'     ## Randomization distribution
#'     p.value <- statmod::permp(A, nperm = n.rep, n1 = m, n2 = n, twosided = (alternative == "two.sided"), method = "auto")
#'   } else if (!randomization) {
#'     ## Permutation distribution
#'     p.value <- A / choose(m + n, m)
#'   }
#'
#'   return(p.value)
#' }
#'
#'
#'
#' ## ----------------------------------------------------------------------------
#' ## Asymmetrically trimmed tests
#' ## ----------------------------------------------------------------------------
#'
#' #' Two-sample location tests based on asymmetrically trimmed means
#' #'
#' #' @description
#' #' \code{asym_trimmed_test()} performs a two-sample location test by using
#' #' asymmetrically trimmed means based on different skewness-
#' #' selector statistics for both samples.
#' #'
#' #'
#' #' @template x
#' #' @template y
#' #' @template type_skewness
#' #' @template alternative
#' #' @template delta
#' #' @template method
#' #' @template n_rep
#' #' @template na_rm
#' #' @template scale_test
#' #'
#' #'
#' #' @details
#' #' The test statistic is an analogue to the one of the ordinary t-test, where
#' #' the difference of the sample means is replaced by the difference of asymmetrically trimmed means
#' #' and the pooled empirical standard deviation is replaced by a pooled winsorized standard deviation
#' #' which uses the asymmetrically trimmed means.
#' #'
#' #' The number of observations trimmed from the lower and the upper end of the sample
#' #' depends on the skewness. Reed & Stark (2004) suggest three possible skewness-selector statistics.
#' #' The argument \code{type} specifies which one is used.
#' #'
#' #' \describe{
#' #' \item{"Q2"}{This skewness-selector statistic compares the difference between the upper and the lower 5\% of
#' #'             the observations to the difference between the upper and lower 50\%.}
#' #' \item{"SK2"}{This skewness-selector statistic compares the difference between the sample minimum and
#' #'              median to the difference between the median and the maximum.}
#' #' \item{"SK5"}{This skewness-selector statistic compares the difference between the sample minimum and
#' #'              mean to the difference between the mean and the maximum.}
#' #' }
#' #'
#' #'
#' #' @return
#' #' A list with class "\code{htest}" containing the following components:
#' #' \item{statistic}{the value of the test statistic.}
#' #' \item{p.value}{the p-value for the test.}
#' #' \item{estimate}{the asymmetrically trimmed means of \code{x} and \code{y}.}
#' #' \item{null.value}{the specified hypothesized value of the mean difference.}
#' #' \item{alternative}{a character string describing the alternative hypothesis.}
#' #' \item{method}{a character string indicating how the p-value was computed.}
#' #' \item{data.name}{a character string giving the names of the data.}
#' #'
#' #' @importFrom Rdpack reprompt
#' #'
#' #' @references
#' #' \insertRef{ReeSta04robu}{robnptests}
#' #'
#' #' @examples
#' #' set.seed(108)
#' #' x <- rnorm(20); y <- rnorm(20)
#' #' asym_trimmed_test(x, y, type = "Q2", method = "asymptotic")
#' #' asym_trimmed_test(x, y, type = "SK2", method = "asymptotic")
#' #' asym_trimmed_test(x, y, type = "SK5", method = "asymptotic")
#' #'
#' #' \dontrun{
#' #' asym_trimmed_test(x, y, type = "SK5", method = "randomization")
#' #' }
#' #'
#' #' @export
#'
#' asym_trimmed_test <- function(x, y, type = c("Q2", "SK2", "SK5"),
#'                               alternative = c("two.sided", "greater", "less"),
#'                               delta = 0,
#'                               method = c("asymptotic", "permutation", "randomization"),
#'                               n.rep = 10000, na.rm = FALSE,
#'                               scale_test = FALSE) {
#'
#'   if (!na.rm & (any(is.na(x)) | any(is.na(y)))) {
#'     return(NA)
#'   } else if (na.rm & (any(is.na(x)) | any(is.na(y)))) {
#'     x <- as.numeric(stats::na.omit(x))
#'     y <- as.numeric(stats::na.omit(y))
#'   }
#'
#'   ## If necessary: Transformation to test for difference in scale
#'   if (scale_test) {
#'     x <- log(x^2)
#'     y <- log(y^2)
#'     delta <- log(delta^2)
#'   }
#'
#'   alternative <- match.arg(alternative)
#'   method <- match.arg(method)
#'   type <- match.arg(type)
#'
#'   ## Error handling
#'   if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
#'     stop ("'delta' must be a single number.")
#'   }
#'
#'   if (!(method %in% c("asymptotic", "permutation", "randomization"))) {
#'     stop (" 'method' must be one of 'asymptotic', 'permutation' or 'randomization'. ")
#'   }
#'
#'   ## Test statistic
#'   t.stat <- asym_trimmed_t(x, y - delta, type = type, na.rm = na.rm)
#'
#'   statistic <- t.stat$statistic
#'   estimates <- t.stat$estimates
#'   df <- t.stat$df
#'
#'   if (method %in% c("permutation", "randomization")) {
#'     ## Calculate permutation distribution
#'     if (method == "randomization") {
#'       randomization <- TRUE
#'     } else {
#'       randomization <- FALSE
#'     }
#'
#'     distribution <- asym_trimmed_perm_distribution(x = x, y = y - delta, type = type,
#'                                                    randomization = randomization, n.rep = n.rep)
#'
#'     ## p-value
#'     p.value <- calc_perm_p_value(statistic, distribution, m = length(x), n = length(y),
#'                                  randomization = randomization, n.rep = n.rep, alternative = alternative)
#'
#'   } else if (method == "asymptotic") {
#'       p.value <- switch (alternative,
#'                          two.sided = 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE),
#'                          greater = stats::pt(statistic, df = df, lower.tail = FALSE),
#'                          less = stats::pt(statistic, df = df, lower.tail = TRUE)
#'     )
#'   }
#'
#'   ## Assign names to results
#'   names(estimates) <- c("Asymmetrically trimmed mean of x", "Asymmetrically trimmed mean of y")
#'   names(delta) <- "location shift"
#'   names(statistic) <- ifelse(scale_test, "S", "D")
#'
#'   if (method == "randomization") {
#'     method = paste("Randomization test based on the", type, "selector statistic")
#'   } else if (method == "permutation") {
#'     method = paste("Exact permutation test based on the", type, "selector statistic")
#'   } else method = paste("Asymptotic test based on the", type, "selector statistic")
#'
#'
#'   dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
#'
#'   res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
#'               estimate = estimates, null.value = delta, alternative = alternative,
#'               method = method, data.name = dname)
#'
#'   class(res) <- "htest"
#'
#'   return(res)
#' }
