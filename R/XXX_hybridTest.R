### EVERYTHING FOR THE HUBER TEST (c-PART)

#' @title Simultaneous Huber-M-estimates of scale and location
#'
#' @description Calculates M-estimates of location and the joined scale of two samples
#'
#' @template x
#' @template y
#' @template k_mest
#'
#' @return Named list containing the following objects
#'         \item{mu.x}{Location estimate of x}
#'         \item{mu.y}{Location estimate of y}
#'         \item{s}{Scale estimate for the joined sample}
#' @import robustbase
#' @export

huber_2sample <- function(x, y, k) {
  m <- length(x)
  n <- length(y)

  N <- m + n

  beta <- 2 * k^2 * (1 - stats::pnorm(k)) + 2 * stats::pnorm(k) - 1 - sqrt(2/pi) * k * exp(-1/2 * k^2)

  s.old <- stats::mad(x) + stats::mad(y)
  #2 * stats::median(c(abs(x - stats::median(x)), abs(y - stats::median(y))))
  mux.old <- stats::median(x)
  muy.old <- stats::median(y)

  #repeat {
  z.x <- (x - mux.old)/s.old
  z.y <- (y - muy.old)/s.old

  s.new <- sqrt(1/((N - 1) * beta) * (sum(robustbase::Mpsi(z.x, psi = "huber", cc = k)^2) + sum(robustbase::Mpsi(z.y, psi = "huber", cc = k)^2)) * s.old^2)

  mux.new <- mux.old + (1/m * sum(robustbase::Mpsi(z.x, psi = "huber", cc = k)) * s.old)/(1/m * sum(robustbase::Mpsi(z.x, psi = "huber", cc = k, deriv = 1)))

  muy.new <- muy.old + (1/n * sum(robustbase::Mpsi(z.y, psi = "huber", cc = k)) * s.old)/(1/n * sum(robustbase::Mpsi(z.y, psi = "huber", cc = k, deriv = 1)))


  #if (abs(mux.new - mux.old) < 1e-6 & abs(muy.new - muy.old) < 1e-6 & abs(s.new - s.old) < 1e-6) {
  #  break
  #}

  s.old <- s.new
  mux.old <- mux.new
  muy.old <- muy.new
  #}

  return(list(mu.x = mux.new, mu.y = muy.new, s = s.new))
}


## ----------------------------------------------------------------------------
## Huber two-sample test
## ----------------------------------------------------------------------------

#' @title Two-sample location tests based on Huber's M-estimator
#'
#' @description Performs a 2-sample location test based on Huber's M-estimator
#'
#' @template x
#' @template y
#' @template delta
#' @template k_huber
#' @template alternative
#' @template var_test
#' @template na_rm
#'
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom for the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the Huber M-estimates of \code{x} and \code{y}.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of test was performed.}
#' \item{data.name}{a character string giving the names of the data.}
#'
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
#' @importFrom Rdpack reprompt
#'
#' @references \insertRef{WeiHot02robu}{robTests}
#'
#' @export

huber_test <- function(x, y, delta = ifelse(var.test, 1, 0), k = 1.8,
                       alternative = c("two.sided", "greater", "less"),
                       var.test = FALSE, na.rm=FALSE) {

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


  if (var.test) {
    names(estimates) <- c("Huber-M of log(x^2)", "Huber-M of log(y^2)")
    names(delta) <- "ratio of variances"
    delta <- exp(delta)
  } else {
    names(estimates) <- c("Huber-M of x", "Huber-M of y")
    names(delta) <- "location shift"
  }

  names(statistic) <- ifelse(var.test, "S", "D")
  method <- "Huber-Test"
  parameter <- m + n - 2
  names(parameter) <- "df"

  res <- list(statistic = statistic, parameter = parameter, p.value = p.value,
              estimate = estimates, null.value = delta,
              alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}


### EVERYTHING FOR THE TRIMMED t-TEST

#' @title Winsorized variance
#'
#' @description
#' \code{win_var} calculates the winsorized variance of a sample.
#'
#'
#' @template x
#' @template gamma_winsorized_variance
#' @template na_rm
#'
#'
#' @return A list containing the following items:
#' \item{var}{The winsorized variance, and}
#' \item{h}{the degrees of freedom used for tests based on trimmed means and the
#' winsorized variance.}
#'
#'
#' @examples
#' ## Generate random sample
#' set.seed(108)
#' x <- rnorm(10)
#'
#' ## Compute 20% winsorized variance
#' win_var(x, gamma = 0.2)
#'
#' @export

win_var <- function(x, gamma = 0, na.rm = FALSE) {
  ## Error handling
  if (gamma < 0 || gamma > 0.5) {
    stop ("gamma has to be in [0, 0.5]")
  }

  ## NA handling
  if (!na.rm & any(is.na(x))) {
    return(NA)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

  ## Calculate winsorized variance
  n <- length(x)

  r <- floor(gamma * n)

  x.sort <- sort(x)
  x.lower <- x.sort[r + 1]
  x.upper <- x.sort[n - r]
  x.sort[x.sort < x.lower] <- x.lower
  x.sort[x.sort > x.upper] <- x.upper

  res <- 1 / (n - 1) * sum((x.sort - mean(x.sort)) ^ 2)
  h <- n - 2 * r

  return(list(var = res, h = h))
}

#' @title Trimmed mean
#'
#' @description
#' \code{trim_mean} calculates a trimmed mean of a sample.
#'
#' @template x
#' @template gamma_trimmed_mean
#' @template na_rm
#'
#' @details
#' This is a wrapper function for the function \code{\link[base]{mean}}.
#'
#' @return
#' The trimmed mean.
#'
#' @examples
#' ## Generate random sample
#' set.seed(108)
#' x <- rnorm(10)
#'
#' ## Compute 20% trimmed mean
#' trim_mean(x, gamma = 0.2)
#'
#' @export

trim_mean <- function(x, gamma = 0.2, na.rm = FALSE) {
  ## Error handling
  if (gamma < 0 || gamma > 0.5) {
    stop ("gamma has to be in [0, 0.5]")
  }

  ## NA handling
  if (!na.rm & any(is.na(x))) {
    return(NA)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

  ## Calculate trimmed mean
  return(mean(x, trim = gamma))
}


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






### TESTS THAT ARE INCLUDED IN THE HYBRID TESTS:

## ----------------------------------------------------------------------------
## Minimum-C-Test
## ----------------------------------------------------------------------------

#' @title Hybrid permutation test for location difference based on the t- and Huber-test
#'
#' @description \code{min_c_test()} performs a hybrid test based on the t-
#'              and Huber-statistic.
#'
#' @template x
#' @template y
#' @template alternative
#' @template delta
#' @template k_minc
#' @template n_rep_hybrid
#' @template na_rm
#' @template n_rep_hybrid
#' @template var_test
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
#' @details The test is introduced in Weichert & Hothorn (2002) and uses the minumum p-value of the t- and the Huber-test as a test statistic.
#' The permutation distribution of the minimum p-value is achieved using the permutation principle according to Efron & Tibshirani (1998).
#'
#' @examples
#' x <- rnorm(10); y <- rnorm(10)
#' min_c_test(x, y, k = 1.8)
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

min_c_test <- function(x, y, alternative = c("two.sided", "greater", "less"), delta = 0,
                       k = 1.8,
                       na.rm = FALSE, n.rep = 1000, var.test = FALSE) {

  alternative <- match.arg(alternative)

  ## If necessary: Transformation to test for difference in scale
  if (var.test) {
    x <- log(x^2)
    y <- log(y^2)
    delta <- log(delta^2)
  }

  z <- c(x, y)
  m <- length(x)
  n <- length(y)

  ## Observed test statistics theta_dach
  t.stat <- stats::t.test(x = x, y = y, mu = delta, alternative = alternative, paired = FALSE, var.equal = TRUE)$statistic
  huber.stat <- huber_test(x = x, y = y, delta = delta, alternative = alternative, k = k)$statistic


  ## Permutation samples and according test statistics theta_dach_stern
  perm.samples <- replicate(n.rep, sample(z))

  perm.t <- apply(perm.samples, 2, function(z) stats::t.test(x = z[1:m], y = z[(m + 1):(m + n)], mu = delta, alternative = alternative, var.equal = TRUE)$statistic)
  perm.huber <- apply(perm.samples, 2, function(z) huber_test(z[1:m], z[(m+1):(m+n)], delta = delta, alternative = alternative, k = k)$statistic)

  pt.obs <- mean(abs(perm.t) >= abs(t.stat))
  phuber.obs <- mean(abs(perm.huber) >= abs(huber.stat))

  obs <- min(pt.obs, phuber.obs)

  statistic <- obs

  pb.t <- sapply(perm.t, function(x) mean(abs(perm.t) >= abs(x)))
  pb.huber <- sapply(perm.huber, function(x) mean(abs(perm.huber) >= abs(x)))

  min.b <- apply(cbind(pb.t, pb.huber), 1, min)

  ## Sachen für die Ausgabe:
  p.value <- mean(min.b <= obs)

  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  names(delta) <- "location shift"
  names(statistic) <- "p"
  method <- "Minimum permutation test based on t- and Huber-statistic"

  res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
              null.value = delta, alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}

## ----------------------------------------------------------------------------
## Minimum-TC-Test
## ----------------------------------------------------------------------------

#' @title Hybrid permutation test for location difference based on the t-, trimmed-t- and Huber-test
#'
#' @description \code{min_tc_test()} performs a hybrid test based on the t-
#'              and Huber-statistic. The significance level is achieved by permutation.
#'
#' @template x
#' @template y
#' @template alternative
#' @template delta
#' @template k_minc
#' @template n_rep_hybrid
#' @template na_rm
#' @template n_rep_hybrid
#' @template var_test
#'
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the minimum p-value of t- and Huber-test.}
#' \item{p.value}{the p-value for the test.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of trimmed t-test was performed.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @details The test is introduced in Weichert & Hothorn (2002) and uses the minumum p-value of the t-, the 20%-
#' trimmed t- and the Huber-test as a test statistic.
#' The permutation distribution of the minimum p-value is achieved using the permutation principle according to Efron & Tibshirani (1998).
#'
#' @examples
#' x <- rnorm(10); y <- rnorm(10)
#' min_tc_test(x, y)
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
#'  \code{\link[robTests]{trimmed_test}}
#'  \code{\link[robTests]{huber_test}}
#' @export
#' @importFrom stats t.test

min_tc_test <- function(x, y, alternative = c("two.sided", "greater", "less"), delta = 0,
                        k = 1.8,
                        na.rm = FALSE, n.rep = 1000, var.test = FALSE) {

  alternative <- match.arg(alternative)

  ## If necessary: Transformation to test for difference in scale
  if (var.test) {
    x <- log(x^2)
    y <- log(y^2)
    delta <- log(delta^2)
  }

  z <- c(x, y)
  m <- length(x)
  n <- length(y)

  ## Observed test statistics theta_dach
  t.stat <- stats::t.test(x = x, y = y, mu = delta, alternative = alternative, paired = FALSE, var.equal = TRUE)$statistic
  t20.stat <- trimmed_test(x = x, y = y, delta = delta, gamma = 0.2, alternative = alternative)$statistic
  huber.stat <- huber_test(x = x, y = y, delta = delta, alternative = alternative, k = k)$statistic

  ## Permutation samples and according test statistics theta_dach_stern
  perm.samples <- replicate(n.rep, sample(z))

  perm.t <- apply(perm.samples, 2, function(z) stats::t.test(x = z[1:m], y = z[(m + 1):(m + n)], mu = delta, alternative = alternative, var.equal = TRUE)$statistic)
  perm.t20 <- apply(perm.samples, 2, function(z) trimmed_test(x = z[1:m], y = z[(m + 1):(m + n)], delta = delta, gamma = 0.2, alternative = alternative)$statistic)
  perm.huber <- apply(perm.samples, 2, function(z) huber_test(z[1:m], z[(m+1):(m+n)], delta = delta, alternative = alternative, k = k)$statistic)

  pt.obs <- mean(abs(perm.t) >= abs(t.stat))
  pt20.obs <- mean(abs(perm.t20) >= abs(t20.stat))
  phuber.obs <- mean(abs(perm.huber) >= abs(huber.stat))

  obs <- min(pt.obs, phuber.obs, pt20.obs)
  statistic <- obs

  pb.t <- sapply(perm.t, function(x) mean(abs(perm.t) >= abs(x)))
  pb.t20 <- sapply(perm.t20, function(x) mean(abs(perm.t20) >= abs(x)))
  pb.huber <- sapply(perm.huber, function(x) mean(abs(perm.huber) >= abs(x)))

  min.b <- apply(cbind(pb.t, pb.huber, pb.t20), 1, min)

  ## Sachen für die Ausgabe:
  p.value <- mean(min.b <= obs)

  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  names(delta) <- "location shift"
  names(statistic) <- "p"
  method <- "Minimum permutation test based on (trimmed) t- and Huber-statistic"

  res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
              null.value = delta, alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}

## ----------------------------------------------------------------------------
## Minimum-t-Test
## ----------------------------------------------------------------------------

#' @title Hybrid permutation test for location difference based on the t- and trimmed t-test
#'
#' @description \code{min_t_test()} performs a hybrid test based on the t-
#'              and Huber-statistic. The significance level is achieved by permutation.
#'
#' @template x
#' @template y
#' @template alternative
#' @template delta
#' @template k_minc
#' @template n_rep_hybrid
#' @template na_rm
#' @template n_rep_hybrid
#' @template var_test
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
#' The permutation distribution of the minimum p-value is achieved using the permutation principle according to Smyth & Phipson (2010).
#'
#' @examples
#' x <- rnorm(10); y <- rnorm(10)
#' min_t_test(x, y)
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{SmyPhi10perm}{robTests}
#' \insertRef{WeiHot02robu}{robTests}
#'
#' @seealso
#'  \code{\link[stats]{t.test}}
#'  \code{\link[robTests]{trimmed_test}}
#'
#' @export

min_t_test <- function(x, y, alternative = c("two.sided", "greater", "less"), delta = 0,
                       na.rm = FALSE, n.rep = 1000, var.test = FALSE) {

  alternative <- match.arg(alternative)

  ## If necessary: Transformation to test for difference in scale
  if (var.test) {
    x <- log(x^2)
    y <- log(y^2)
    delta <- log(delta^2)
  }

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





#' @title Hybrid permutation test for location differences
#'
#' @description \code{hybrid_test()} performs a hybrid test based on the p-values of different test statistics
#'
#' @template x
#' @template y
#' @template type_hybrid
#' @template alternative
#' @template delta
#' @template k_hybrid
#' @template na_rm
#' @template n_rep_hybrid
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
