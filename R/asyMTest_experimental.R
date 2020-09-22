## ----------------------------------------------------------------------------
## Two-sample location test based on difference of M-estimators
## ----------------------------------------------------------------------------

#' Two-sample location tests based on difference of M-estimators
#'
#' @description
#' \code{m_test} performs a two-sample location test based on
#' the difference of M-estimates of both samples.
#'
#' @template x
#' @template y
#' @template alternative
#' @template delta
#' @template method
#' @template psi
#' @template k_mest
#' @template n_rep
#' @template na_rm
#' @template var_test
#' @template wobble
#' @template seed
#'
#' @export

m_test <- function(x, y, alternative = c("two.sided", "greater", "less"),
                   delta = ifelse(var.test, 1, 0),
                   method = c("asymptotic", "permutation", "randomization"),
                   psi = c("huber", "hampel", "bisquare"),
                   k = .Mpsi.tuning.default(psi),
                   n.rep = 10000, na.rm = FALSE,
                   var.test = FALSE, wobble = FALSE, wobble.seed = NULL) {

  alternative <- match.arg(alternative)
  #method <- match.arg(method)
  psi <- match.arg(psi)

  ## Names of data sets
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  ## NA handling
  if (!na.rm & (any(is.na(x)) | any(is.na(y)))) {
    return(NA)
  } else if (na.rm & (any(is.na(x)) | any(is.na(y)))) {
    x <- as.numeric(stats::na.omit(x))
    y <- as.numeric(stats::na.omit(y))
  }

  ## Check sample sizes
  if (length(x) < 5 || length(y) < 5) {
    stop("Both samples need at least 5 non-missing values.")
  }

  ## Add random noise if user wants to wobble
  if (wobble) {

    if (is.null(wobble.seed)) wobble.seed <- sample(1e6, 1)
    set.seed(wobble.seed)

    xy <- wobble(x, y)
    x <- xy$x
    y <- xy$y

    warning(paste0("Added random noise to x and y. The seed is ",
                   wobble.seed, "."))
  }

  ## If necessary: Transformation to test for difference in scale
  if (var.test) {
    if (any(c(x, y) == 0)) {

      if (is.null(wobble.seed)) wobble.seed <- sample(1e6, 1)
      set.seed(wobble.seed)

      xy <- wobble(x, y, check = FALSE)
      x <- xy$x
      y <- xy$y

      warning(paste0("Added random noise before log transformation due to zeros in the sample. The seed is ",
                     wobble.seed, "."))
    }
    x <- log(x^2)
    y <- log(y^2)
    delta <- log(delta^2)
  }

  ## Error handling
  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
    stop ("'delta' must be a single number.")
  }

  if (!all((method %in% c("asymptotic", "permutation", "randomization")))) {
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
    ## _________________________________________________________________________
    ## Compute test statistic
    ## _________________________________________________________________________
    stats <- m_test_statistic(x, y - delta, psi = psi, k = k)
    statistic <- stats$statistic

    estimates <- stats$estimates

    distribution <- mest_perm_distribution(x = x, y = y - delta, randomization = (method == "randomization"),
                                           n.rep = n.rep, psi = psi, k = k)

    p.value <- calc_perm_p_value(statistic, distribution, m = length(x), n = length(y),
                                 randomization = randomization, n.rep = n.rep, alternative = alternative)
  } else if (method == "asymptotic") {
    m <- length(x)
    n <- length(y)

    ## _________________________________________________________________________
    ## Compute test statistic
    ## _________________________________________________________________________

    ## M-estimates for both samples
    est.x <- m_est(x = x, psi = psi, k = k)$est
    est.y <- m_est(x = y, psi = psi, k = k)$est

    ## Estimator for \nu
    psi.x <- robustbase::Mpsi((x - est.x)/mad(x), psi = psi, cc = k)
    rho.x <- robustbase::Mpsi((x - est.x)/mad(x), psi = psi, cc = k, deriv = 1)

    psi.y <- robustbase::Mpsi((y - est.y)/mad(y), psi = psi, cc = k)
    rho.y <- robustbase::Mpsi((y - est.y)/mad(y), psi = psi, cc = k, deriv = 1)

    nu.x <- mean(psi.x^2)/(mean(rho.x)^2)
    nu.y <- mean(psi.y^2)/(mean(rho.y)^2)

    ## Value of the test statistic
    statistic <- (est.x - est.y) /
      sqrt((n * robustbase::scaleTau2(x, consistency = FALSE)^2 * nu.x +
              m * robustbase::scaleTau2(y, consistency = FALSE)^2 * nu.y) / (m * n))

    # ## ___________________________________________________________________________
    # ## Compute test statistic
    # ## ___________________________________________________________________________

    # ## M-estimates
    # est.x <- robustbase::huberM(x, k = 1.5)
    # est.y <- robustbase::huberM(y, k = 1.5)
    #
    # ## Estimator for \nu
    # psi.x <- robustbase::Mpsi((x - est.x$mu)/est.x$s, psi = "huber", cc = 1.5)
    # rho.x <- robustbase::Mpsi((x - est.x$mu)/est.x$s, psi = "huber", cc = 1.5, deriv = 1)
    #
    # psi.y <- robustbase::Mpsi((y - est.y$mu)/est.y$s, psi = "huber", cc = 1.5)
    # rho.y <- robustbase::Mpsi((y - est.y$mu)/est.y$s, psi = "huber", cc = 1.5, deriv = 1)
    #
    # #nu.x <- 2 * (k^2 * (1 - pnorm(k)) + pnorm(k) - 0.5 - k * dnorm(k))/(pnorm(k) - pnorm(-k))^2 #
    # nu.x <- mean(psi.x^2)/(mean(rho.x)^2)
    # nu.y <- mean(psi.y^2)/(mean(rho.y)^2)
    #
    # #res[i] <- sqrt((m * n)/(m + n)) * (est.x$mu - est.y$mu)/((nu.x + nu.y)/2)
    # statistic <- (est.x$mu - est.y$mu) /
    #   sqrt((n * robustbase::scaleTau2(x, consistency = FALSE)^2 * nu.x +
    #           m * robustbase::scaleTau2(y, consistency = FALSE)^2 * nu.y) / (m * n))
    # #
    # ## ___________________________________________________________________________
    # ## Test decision
    # ## ___________________________________________________________________________
    # p.value <- switch (alternative,
    #                    two.sided = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
    #                    greater = stats::pnorm(statistic, lower.tail = FALSE),
    #                    less = stats::pnorm(statistic, lower.tail = TRUE)
    # )
  }

  # if (var.test) {
  #   names(estimates) <- c("M-est. of log(x^2)", "M-est. of log(y^2)")
  #   names(delta) <- "ratio of variances"
  #   delta <- exp(delta)
  # } else {
  #   names(estimates) <- c("M-est. of x", "M-est. of y")
  #   names(delta) <- "location shift"
  # }
  # names(statistic) <- ifelse(var.test, "S", "D")
  #
  # if (method == "randomization") {
  #   method = paste("Randomization test based on the ", paste0(toupper(substring(psi, 1, 1)), substring(psi, 2, nchar(psi))), "M-estimator")
  # } else method = paste("Exact permutation test based on the", psi, "M-estimator")

  res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
              estimate = estimates, null.value = delta, alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}

res <- numeric(1000)

set.seed(108)
for (i in 1:1000) {
  x <- rnorm(500)
  y <- rnorm(500)

  res[i] <- m_test(x = x, y = y, alternative = "two.sided", method = "asymptotic")$statistic

}


res <- numeric(10000)

m <- n <- 100


for (i in 1:10000) {
  x <- rt(m, df = 3)
  y <- rt(n, df = 3)

  ## M-estimates for both samples
  est.x <- m_est(x = x, psi = psi, k = k)
  est.y <- m_est(x = y, psi = psi, k = k)

  ## Estimator for \nu
  psi.x1 <- robustbase::Mpsi((x - est.x$est)/mad(x), psi = psi, cc = k)
  rho.x <- robustbase::Mpsi((x - est.x$est)/mad(x), psi = psi, cc = k, deriv = 1)

  psi.y <- robustbase::Mpsi((y - est.y$est)/mad(y), psi = psi, cc = k)
  rho.y <- robustbase::Mpsi((y - est.y$est)/mad(y), psi = psi, cc = k, deriv = 1)

  nu.x <- mean(psi.x^2)/(mean(rho.x)^2)
  nu.y <- mean(psi.y^2)/(mean(rho.y)^2)

  #res[i] <- sqrt((m * n)/(m + n)) * (est.x$mu - est.y$mu)/((nu.x + nu.y)/2)
  res[i] <- (est.x$est - est.y$est) /
    sqrt((n * robustbase::scaleTau2(x, consistency = FALSE)^2 * nu.x +
            m * robustbase::scaleTau2(y, consistency = FALSE)^2 * nu.y) / (m * n))
}




res <- numeric(10000)

m <- n <- 100


for (i in 1:10000) {
  x <- rnorm(m)
  y <- rnorm(n)

  ## M-estimates
  est.x <- huberM(x, k = 1.5)
  est.y <- huberM(y, k = 1.5)

  ## Estimator for \nu
  psi.x <- robustbase::Mpsi((x - est.x$mu)/est.x$s, psi = "huber", cc = 1.5)
  rho.x <- robustbase::Mpsi((x - est.x$mu)/est.x$s, psi = "huber", cc = 1.5, deriv = 1)

  psi.y <- robustbase::Mpsi((y - est.y$mu)/est.y$s, psi = "huber", cc = 1.5)
  rho.y <- robustbase::Mpsi((y - est.y$mu)/est.y$s, psi = "huber", cc = 1.5, deriv = 1)

  #nu.x <- 2 * (k^2 * (1 - pnorm(k)) + pnorm(k) - 0.5 - k * dnorm(k))/(pnorm(k) - pnorm(-k))^2 #
  nu.x <- mean(psi.x^2)/(mean(rho.x)^2)
  nu.y <- mean(psi.y^2)/(mean(rho.y)^2)

  #res[i] <- sqrt((m * n)/(m + n)) * (est.x$mu - est.y$mu)/((nu.x + nu.y)/2)
  res[i] <- (est.x$mu - est.y$mu)/sqrt((n * scaleTau2(x, sigma0=mad(x))^2 * nu.x + m * scaleTau2(y, sigma0 = mad(y))^2 * nu.y)/(m * n))
}


hist(res, freq = FALSE, xlim = c(-5, 5), ylim = c(0, 0.5))
curve(dnorm(x, sd = 1), from = -5, to = 5, add = TRUE)


qqnorm(res)


