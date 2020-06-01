## ----------------------------------------------------------------------------
## Calculation permutation distribution
## ----------------------------------------------------------------------------

#' Permutation Distribution
#'
#' \code{perm_distribution()} calculates the permutation distribution for
#' several test statistics.
#'
#' @template x
#' @template y
#' @template type_rob_perm
#' @template randomization
#' @template n_rep
#'
#' @return Vector with permutation distribution.
#'
#' @details see documentation of \code{rob_perm_statistic()} for a description of
#' the \code{type}-parameter
#'
#' @export

perm_distribution <- function(x, y, type, randomization = FALSE, n.rep = 10000) {
  ## Sample sizes
  m <- length(x)
  n <- length(y)

  ## Error handling
  if (randomization & n.rep > choose(m + n, m)) {
    stop (paste0("'n.rep' must not be larger than ", choose(m + n, m), ", the number of all splits."))
  }

  ## Splits in two samples
  if (!randomization) {
    ## Computation of permutation distribution

    complete <- c(x, y)
    splits <- gtools::combinations((m + n), m, 1:(m + n))

    distribution <- apply(splits, 1, function(s) rob_perm_statistic(x = complete[s], y = complete[-s], type)$statistic)
  } else if (randomization) {
    ## Computation of randomization distribution

    splits <- replicate(n.rep, sample(c(x, y)))

    distribution <- apply(splits, 2, function(s) rob_perm_statistic(x = s[1:m], y = s[(m + 1):(m + n)], type)$statistic)
  }

  return(distribution)
}

#' @title Permutation distribution for M-statistics
#'
#' @description \code{mest_perm_distribution} calculates the permutation distribution for M-statistics from
#'              \code{m_test_statistic}.
#'
#' @template x
#' @template y
#' @template psi
#' @template k1_mest
#' @template randomization
#' @template n_rep_m_test
#'
#' @return Vector with permutation distribution.
#'
#' @export

mest_perm_distribution <- function(x, y, psi, k1, randomization = FALSE, n.rep = NULL) {
  m <- length(x)
  n <- length(y)

  if (!randomization) {
    complete <- c(x, y)
    splits <- gtools::combinations((m + n), h, 1:(m + n))

    distribution <- apply(splits, 1, function(s) m_test_statistic(x = complete[s], y = complete[-s],
                                                                  psi = psi, k = k1)$statistic)
  }
  else if (randomization) {
    splits <- replicate(n.rep, sample(c(x, y)))

    distribution <- apply(splits, 2, function(s) m_test_statistic(x = s[1:m], y = s[(m + 1):(m + n)], psi = psi, k = k1)$statistic)
  }

  return(distribution)
}

#' @title Permutation distribution for asymmetrically trimmed statistics
#'
#' @description \code{asym_trimmed_perm_distribution} calculates the permutation distribution for the asymmetrically
#' trimmed statistics from \code{asym_trimmed_test}.
#'
#' @template x
#' @template y
#' @template type_skewness
#' @template randomization
#' @template n_rep
#'
#' @return Vector with permutation distribution.
#'
#' @export

asym_trimmed_perm_distribution <- function(x, y, type, randomization = FALSE, n.rep = NULL) {
  ## Sample sizes
  m <- length(x)
  n <- length(y)

  ## Splits in two samples
  if (!randomization) {
    complete <- c(x, y)
    splits <- gtools::combinations((m + n), m, 1:(m + n))

    distribution <- apply(splits, 1, function(s) asym_trimmed_t(x = complete[s], y = complete[-s], type)$statistic)

  } else if (randomization) {
    splits <- replicate(n.rep, sample(c(x, y)))

    distribution <- apply(splits, 2, function(s) asym_trimmed_t(x = s[1:m], y = s[(m + 1):(m + n)], type)$statistic)
  }

  return(distribution)
}


## ----------------------------------------------------------------------------
## Calculate p-value for permutation tests
## ----------------------------------------------------------------------------

#' Calculation of permutation p-value
#'
#' @description
#' \code{calc_perm_p_value} calculates the permutation p-value following Phipson & Smyth (2010).
#'
#' @template statistic
#' @template distribution
#' @template m
#' @template n
#' @template randomization
#' @template n_rep
#' @template alternative
#'
#' @return
#' p.value for the specified alternative.
#'
#' @references
#' \insertRef{SmyPhi10perm}{robTests}
#'
#' @export

calc_perm_p_value <- function(statistic, distribution, m, n, randomization, n.rep, alternative) {

  ## Number of permutations leading to test statistic at least as extreme
  ## as observed
  A <- switch(alternative,
              two.sided = sum(abs(distribution) >= abs(statistic)),
              greater = sum(distribution >= statistic),
              less = sum(distribution <= statistic)
  )

  ## Computation of p-value
  if (randomization) {
    ## Randomization distribution
    p.value <- statmod::permp(A, nperm = n.rep, n1 = m, n2 = n, twosided = (alternative == "two.sided"), method = "auto")
  } else if (!randomization) {
    ## Permutation distribution
    p.value <- A / choose(m + n, m)
  }

  return(p.value)
}

