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
#' @keywords internal

perm_distribution <- function(x, y, type, randomization = FALSE, n.rep = 10000) {

  ## Check input arguments
  stopifnot("'x' is missing." = !missing(x))
  stopifnot("'y' is missing." = !missing(y))
  stopifnot("'type' is missing." = !missing(type))

  checkmate::assert_numeric(x, min.len = 5, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(y, min.len = 5, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_choice(type, choices = c("S1", "S2", "S3", "S4"), null.ok = FALSE)
  checkmate::assert_flag(randomization, na.ok = FALSE, null.ok = FALSE)
  checkmate::assert_count(n.rep, na.ok = FALSE, positive = TRUE, null.ok = FALSE)

  ## Sample sizes
  m <- length(x)
  n <- length(y)

  ## For the randomization distribution, the value of 'n.rep' is bounded by the
  ## number of possible splits into two samples
  if (randomization & (n.rep > choose(m + n, m))) {
    stop (paste0("'n.rep' must not be larger than ", choose(m + n, m), ", the number of all splits."))
  }

  ## Splits in two samples
  if (!randomization) {
    ## Computation of the permutation distribution

    complete <- c(x, y)
    splits <- gtools::combinations((m + n), m, 1:(m + n))

    distribution <- apply(splits, 1, function(s) rob_perm_statistic(x = complete[s], y = complete[-s], type)$statistic)
  } else if (randomization) {
    ## Computation of the randomization distribution

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
#' @keywords internal

mest_perm_distribution <- function(x, y, psi, k, randomization = FALSE, n.rep = 10000) {

  ## Check input arguments
  stopifnot("'x' is missing." = !missing(x))
  stopifnot("'y' is missing." = !missing(y))
  stopifnot("'psi' is missing." = !missing(k))
  stopifnot("'k' is missing." = !missing(psi))

  checkmate::assert_numeric(x, min.len = 5, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(y, min.len = 5, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_choice(psi, choices = c("huber", "hampel", "tukey"), null.ok = FALSE)
  checkmate::assert_number(k, na.ok = FALSE, lower = 0, finite = TRUE, null.ok = FALSE)
  checkmate::assert_flag(randomization, na.ok = FALSE, null.ok = FALSE)
  checkmate::assert_count(n.rep, na.ok = FALSE, positive = TRUE, null.ok = FALSE)

  ## Sample sizes
  m <- length(x)
  n <- length(y)

  ## For the randomization distribution, the value of 'n.rep' is bounded by the
  ## number of possible splits into two samples
  if (randomization & (n.rep > choose(m + n, m))) {
    stop (paste0("'n.rep' must not be larger than ", choose(m + n, m), ", the number of all splits."))
  }

  ## Splits in two samples
  if (!randomization) {
    ## Computation of the permutation distribution
    complete <- c(x, y)
    splits <- gtools::combinations((m + n), m, 1:(m + n))

    distribution <- apply(splits, 1, function(s) m_test_statistic(x = complete[s], y = complete[-s],
                                                                  psi = psi, k = k)$statistic)
  } else if (randomization) {
    ## Computation of the randomization distribution
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
#' @keywords internal

calc_perm_p_value <- function(statistic, distribution, m, n, randomization, n.rep, alternative) {

  ## Check input arguments
  stopifnot("'statistic' is missing." = !missing(statistic))
  stopifnot("'distribution' is missing." = !missing(distribution))
  stopifnot("'m' is missing." = !missing(m))
  stopifnot("'n' is missing." = !missing(n))
  stopifnot("'randomization' is missing." = !missing(randomization))
  stopifnot("'n.rep' is missing." = !missing(n.rep))
  stopifnot("'alternative' is missing." = !missing(alternative))

  checkmate::assert_number(statistic, na.ok = FALSE, finite = TRUE, null.ok = FALSE)
  checkmate::assert_numeric(distribution, finite = TRUE, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_count(m, na.ok = FALSE, positive = TRUE, null.ok = FALSE)
  checkmate::assert_count(n, na.ok = FALSE, positive = TRUE, null.ok = FALSE)
  checkmate::assert_flag(randomization, na.ok = FALSE, null.ok = FALSE)
  checkmate::assert_count(n.rep, na.ok = FALSE, positive = TRUE, null.ok = FALSE)
  checkmate::assert_choice(alternative, choices = c("two.sided", "greater", "less"), null.ok = FALSE)

  ## Number of permutations leading to test statistic at least as extreme
  ## as the observed value
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
