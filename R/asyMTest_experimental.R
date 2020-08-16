m_test2 <- function(x, y, alternative = c("two.sided", "greater", "less")) {

  alternative <- match.arg(alternative)

  m <- length(x)
  n <- length(y)

  ## ___________________________________________________________________________
  ## Compute test statistic
  ## ___________________________________________________________________________
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
  statistic <- (est.x$mu - est.y$mu)/sqrt((n * scaleTau2(x, consistency = FALSE)^2 * nu.x + m * scaleTau2(y, consistency = FALSE)^2 * nu.y)/(m * n))

  ## ___________________________________________________________________________
  ## Test decision
  ## ___________________________________________________________________________
  p.value <- switch (alternative,
                     two.sided = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
                     greater = stats::pnorm(statistic, lower.tail = FALSE),
                     less = stats::pnorm(statistic, lower.tail = TRUE)
  )

  return(list(statistic = statistic, p.value = p.value))

}

res <- numeric(1000)

set.seed(108)
for (i in 1:1000) {
  x <- rt(500, df = 3)
  y <- rt(500, df = 3)

  res[i] <- m_test2(x = x, y = y, alternative = "two.sided")$p.value

}



res <- numeric(1000)

m <- n <- 100


for (i in 1:1000) {
  x <- rt(m, df = 3)
  y <- rt(n, df = 3)

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
  res[i] <- (est.x$mu - est.y$mu)/sqrt((n * scaleTau2(x)^2 * nu.x + m * scaleTau2(y)^2 * nu.y)/(m * n))
}


hist(res, freq = FALSE, xlim = c(-5, 5), ylim = c(0, 0.4))
curve(dnorm(x, sd = 1), from = -5, to = 5, add = TRUE)

qqnorm(res)


