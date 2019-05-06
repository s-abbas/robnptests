## ----------------------------------------------------------------------------
## Calculation permutation distribution
## ----------------------------------------------------------------------------

#' Permutation Distribution
#'
#' \code{perm_distribution()} calculates the permutation distribution for
#' several test statistics.
#'
#' @inheritParams rob_perm_statistic
#' @param sampled logical value indicating whether all splits should be
#' considered or only a random sample.
#' @param n.rep integer specifying the sample size in case of taking a random
#' sample to approximate the permutation distribution.
#'
#' @return Vector with permutation distribution.
#'
#' @details see documentation of \code{rob_perm_statistic()} for a description of
#' the \code{type}-parameter
#'
#' @export

perm_distribution <- function(x, y, type, sampled = FALSE, n.rep = NULL) {
  ## Sample sizes
  h <- length(x)
  k <- length(y)

  ## Splits in two samples
  if (!sampled) {
    complete <- c(x, y)
    splits <- gtools::combinations((h + k), h, 1:(h + k))

    distribution <- apply(splits, 1, function(s) rob_perm_statistic(x = complete[s], y = complete[-s], type)$statistic)

  } else if (sampled) {
    splits <- replicate(n.rep, sample(c(x, y)))

    distribution <- apply(splits, 2, function(s) rob_perm_statistic(x = s[1:h], y = s[(h + 1):(h + k)], type)$statistic)
  }

  return(distribution)
}

#' @title Permutation distribution for M-statistics
#'
#' @description \code{mest_perm_distribution} calculates the permutation distribution for M-statistics from
#'              \code{m_test_statistic}.
#'
#' @inheritParams m_test_statistic
#' @param psi kernel used for optimization.
#' @param k1 tuning parameter(s) for the respective kernel function.
#' @param sampled logical value indicating if the exact permutation distribution (\code{FALSE})
#'                or only \code{n.rep} samples drawn with replacement (\code{TRUE})
#'                should be used. Default is FALSE.
#' @param n.rep integer specifying the sample size in case of taking a random
#'              sample to approximate the permutation distribution.
#'
#' @return Vector with permutation distribution.
#'
#' @export

mest_perm_distribution <- function(x, y, psi, k1, sampled = FALSE, n.rep = NULL) {
  h <- length(x)
  k <- length(y)

  if (!sampled) {
    complete <- c(x, y)
    splits <- gtools::combinations((h + k), h, 1:(h + k))

    distribution <- apply(splits, 1, function(s) m_test_statistic(x = complete[s], y = complete[-s],
                                                                  psi = psi, k = k1)$statistic)
  }
  else if (sampled) {
    splits <- replicate(n.rep, sample(c(x, y)))

    distribution <- apply(splits, 2, function(s) m_test_statistic(x = s[1:h], y = s[(h+1):(h+k)], psi = psi, k = k1)$statistic)
  }

  return(distribution)
}

#' @title Permutation distribution for asymmetrically trimmed statistics
#'
#' @description \code{asym_trimmed_perm_distribution} calculates the permutation distribution for the asymmetrically
#' trimmed statistics from \code{asym_trimmed_test}.
#'
#' @inheritParams asym_trimmed_test
#'
#' @param sampled logical value indicating if the exact permutation distribution (\code{FALSE})
#'                or only \code{n.rep} samples drawn with replacement (\code{TRUE})
#'                should be used. Default is FALSE.
#' @param n.rep integer specifying the sample size in case of taking a random
#'              sample to approximate the permutation distribution.
#'
#' @return Vector with permutation distribution.
#'
#' @export

asym_trimmed_perm_distribution <- function(x, y, type, sampled = FALSE, n.rep = NULL) {
  ## Sample sizes
  h <- length(x)
  k <- length(y)

  ## Splits in two samples
  if (!sampled) {
    complete <- c(x, y)
    splits <- gtools::combinations((h + k), h, 1:(h + k))

    distribution <- apply(splits, 1, function(s) asym_trimmed_t(x = complete[s], y = complete[-s], type)$statistic)

  } else if (sampled) {
    splits <- replicate(n.rep, sample(c(x, y)))

    distribution <- apply(splits, 2, function(s) asym_trimmed_t(x = s[1:h], y = s[(h + 1):(h + k)], type)$statistic)
  }

  return(distribution)
}


## ----------------------------------------------------------------------------
## Calculate p-value for permutation tests
## ----------------------------------------------------------------------------

#' Calculation of permutation p-value
#'
#' @description
#' \code{calc_perm_p_value} calculates the permutation p-value following Smyth & Phipson (2010).
#'
#' @param statistic observed value of the test statistic.
#' @param distribution numeric vector with permutation distribution.
#' @param m integer value giving size of first sample.
#' @param n integer value giving size of second sample.
#' @param sampled logical value specifying whether the distribution consists of
#'                all permutations or only a random subset drawn with replacement.
#' @param n.rep integer value giving the number of replications
#'              (only needed if \code{sampled = TRUE}).
#' @param alternative character string specifying the alternative hypothesis,
#'                    must be one of "\code{two.sided}" (default),
#'                    "\code{greater}" or "\code{less}".
#'
#' @return
#' p.value for the specified alternative.
#'
#' @references
#' \insertRef{SmyPhi10perm}{robTests}
#'
#' @export

calc_perm_p_value <- function(statistic, distribution, m, n, sampled, n.rep, alternative) {

  ## Number of permutations leading to test statistic at least as extreme
  ## as observed
  A <- switch(alternative,
              two.sided = sum(abs(distribution) >= abs(statistic)),
              greater = sum(distribution >= statistic),
              less = sum(distribution <= statistic)
  )

  if (sampled) {
    ## Random subset of permutations with replacement
    p.value <- statmod::permp(A, nperm = n.rep, n1 = m, n2 = n, twosided = (alternative == "two.sided"))
  } else if (!sampled) {
    ## All permutations
    p.value <- A / choose(m + n, m)
  }

  return(p.value)
}

