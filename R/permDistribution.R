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
#' @return Vector with permutation distribution of the test statistic specified by \code{type}.
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
#' @description \code{mest_perm_distribution} calculates the permutation distribution for the M-statistics from
#'              \code{m_test_statistic}.
#'
#' @template x
#' @template y
#' @template psi
#' @template k_mest
#' @template randomization
#' @template n_rep
#'
#' @return Vector with permutation distribution of the test statistic specified by \code{psi}
#'         and \code{k}.
#'
#' @references
#' \insertRef{MaeRouCro20robu}{robTests}
#'
#' @export

mest_perm_distribution <- function(x, y, psi, k, randomization = FALSE, n.rep = NULL) {
  m <- length(x)
  n <- length(y)

  if (!randomization) {
    complete <- c(x, y)
    splits <- gtools::combinations((m + n), m, 1:(m + n))

    distribution <- apply(splits, 1, function(s) m_test_statistic(x = complete[s], y = complete[-s],
                                                                  psi = psi, k = k)$statistic)
  } else if (randomization) {
    splits <- replicate(n.rep, sample(c(x, y)))

    distribution <- apply(splits, 2, function(s) m_test_statistic(x = s[1:m], y = s[(m + 1):(m + n)], psi = psi, k = k)$statistic)
  }

  return(distribution)
}

## ----------------------------------------------------------------------------
## Calculate p-value for permutation tests
## ----------------------------------------------------------------------------

#' Calculation of permutation p-value
#'
#' @description
#' \code{calc_perm_p_value} calculates the permutation p-value following \insertCite{PhiSmy10perm;textual}{robTests}.
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
#' p-value for the specified alternative.
#'
#' @references
#' \insertRef{PhiSmy10perm}{robTests}
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
