#' @title Preprocess data for the robust two sample tests
#'
#' @description
#' \code{preprocess_data} is a helper function that performs several
#' preprocessing steps on the data before performing the two-sample tests.
#'
#' @template x
#' @template y
#' @template delta
#' @template na_rm
#' @template wobble
#' @template wobble_seed
#' @template var_test
#'
#' @details
#' The preprocessing steps include the removal of missing values and, if
#' specified, wobbling and a transformation of the observations to test for
#' differences in scale.
#'
#' @return A named list containing the following components:
#'         \item{\code{x}}{the (possibly transformed) input vector \code{x}.}
#'         \item{\code{y}}{the (possibly transformed) input vector \code{y}.}
#'         \item{\code{delta}}{the (possibly transformed) input value
#'                             \code{delta}.}
#'
#' @keywords internal

preprocess_data <- function(x, y, delta, na.rm, wobble, wobble.seed, var.test) {

  ## Remove missing values ----
  if (na.rm) {
    x <- as.vector(na.omit(x))
    y <- as.vector(na.omit(y))
  }

  ## Check sample sizes ----
  if (length(x) < 5 || length(y) < 5) {
    stop("Both samples need at least 5 non-missing values.")
  }

  ## Wobbling ----
  if (wobble) {
    # Random noise is only added, if 'wobbling' = TRUE and there is at least
    # one duplicated value in the single samples or the joint sample
    if (!(length(unique(x)) == length(x) &
          length(unique(y)) == length(y) &
          length(unique(c(x, y))) == length(c(x, y)))) {

      # Set seed for generating random noise
      if (missing(wobble.seed) | is.null(wobble.seed)) {
        wobble.seed <- sample(1e6, 1)
      }
      set.seed(wobble.seed)

      # Add random noise
      xy <- wobble(x, y, check = FALSE)
      x <- xy$x
      y <- xy$y

      warning(paste0("Added random noise to x and y. The seed is ",
                     wobble.seed, "."))
    }
  }

  ## Transformation of observations for variance test ----
  if (var.test) {

    if (delta == 0) {
      ## 'delta' needs to be larger than zero
      stop("The logarithm of 0 is not defined. Please use another value for 'delta'.")
    }

    if (any(c(x, y) == 0)) {
      # The log-transformation only works for non-zero values
      if (missing(wobble.seed) | is.null(wobble.seed)) wobble.seed <- sample(1e6, 1)

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
    ## The transformation of Delta needs to be done in the function!!
  }

  return(list(x = x, y = y, delta = delta))
}

#' @title Checks for input arguments
#'
#' @description
#' \code{check_test_input} is a helper functions that contains checks for the
#' input arguments of the two-sample tests.
#'
#' @template x
#' @template y
#' @template alternative
#' @template delta
#' @param scale a character string specifying the scale estimator used for
#'              standardization in the test statistic; must be one of \code{"S1"},
#'              \code{"S2"}, \code{"S3"}, and \code{"S4"}.
#' @template n_rep
#' @template na_rm
#' @template var_test
#' @template wobble
#' @template wobble_seed
#' @template gamma_trimmed_test
#' @template psi
#' @template k_mest
#' @template test_name
#'
#' @details
#' The two-sample tests in this package share similar arguments. To reduce the
#' amount of repetitive code, this function contains the argument checks so that
#' only \code{check_test_input} needs to be called within the functions for
#' the two-sample tests.
#'
#' The scale estimators \code{"S1"} and \code{"S2"} can only be used in
#' combination with \code{test.name = "hl1_test"} or \code{test.name = "hl2_test"}.
#' The estimators \code{"S3"} and \code{"S4"} can only be used with
#' \code{test.name = "med_test"}.
#'
#' @return An error message if a check fails.
#'
#' @keywords internal

check_test_input <- function(x,
                             y,
                             alternative,
                             delta,
                             method,
                             scale,
                             n.rep,
                             na.rm,
                             var.test,
                             wobble,
                             wobble.seed,
                             gamma = NULL,
                             psi = NULL,
                             k = NULL,
                             test.name) {

  ## Checks that are necessary for all tests ----
  checkmate::assert_numeric(x, finite = TRUE, all.missing = FALSE, min.len = 5, null.ok = FALSE)
  checkmate::assert_numeric(y, finite = TRUE, all.missing = FALSE, min.len = 5, null.ok = FALSE)
  checkmate::assert_subset(alternative, choices = c("two.sided", "greater", "less"), empty.ok = FALSE)
  checkmate::assert_subset(method, choices = c("asymptotic", "permutation", "randomization"), empty.ok = FALSE)
  checkmate::assert_count(n.rep, na.ok = FALSE, positive = TRUE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)
  checkmate::assert_flag(var.test, na.ok = FALSE, null.ok = FALSE)
  if (var.test) {
    checkmate::assert_numeric(delta, lower = 1e-6, finite = TRUE, any.missing = FALSE, len = 1, null.ok = FALSE)
  } else if (!var.test) {
    checkmate::assert_numeric(delta, finite = TRUE, any.missing = FALSE, len = 1, null.ok = FALSE)
  }
  checkmate::assert_numeric(wobble.seed, finite = TRUE, any.missing = FALSE, len = 1, null.ok = TRUE)

  ## Checks for HL1-, HL2-, and MED-tests ----
  if (test.name %in% c("hl1_test", "hl2_test")) {
    checkmate::assert_subset(scale, choices = c("S1", "S2"), empty.ok = FALSE)
  } else if (test.name == "med_test") {
    checkmate::assert_subset(scale, choices = c("S3", "S4"), empty.ok = FALSE)
  }
  if (!(test.name %in% c("trimmed_test", "m_test"))) {
    checkmate::assert_flag(wobble, na.ok = FALSE, null.ok = FALSE)
  }

  ## Additional checks for trimmed t-test ---
  if (test.name == "trimmed_test") {
    checkmate::assert_number(gamma, na.ok = FALSE, lower = 0, upper = 0.5, finite = TRUE, null.ok = FALSE)
  }

  ## Additional checks for M-tests ----
  if (test.name == "m_test") {
    checkmate::assert_subset(psi, choices = c("huber", "hampel", "bisquare"), empty.ok = FALSE)
    checkmate::assert_numeric(k, lower = 0, len = ifelse(psi == "hampel", 3, 1), finite = TRUE, any.missing = FALSE, null.ok = FALSE)
  }
}

# #' @title Remove missing values in \code{x} and \code{y}
# #'
# #' @description
# #' \code{remove_missing_values} is a helper function that removes missing values
# #' from the input sample vectors.
# #'
# #' @template x
# #' @template y
# #' @template na_rm
# #'
# #' @details
# #' Both samples need at least five non-NA observations.
# #'
# #' @return
# #' If one of the sample contains a missing value and \code{na.rm = FALSE},
# #' \code{NA_real_} is returned. If the number of non-missing observations
# #' in one of the samples is less than five, the function stops with an error
# #' message. Otherwise, the function returns a named list
# #' containing the following compoents:
# #' \item{\code{x}}{input vector \code{x} without missing values.}
# #' \item{\code{y}}{input vector \code{y} without missing values.}
# #'
# #' @keywords internal
#
# remove_missing_values <- function(x, y) {
#
#   # Remove missing values ----
#   x <- as.vector(stats::na.omit(x))
#   y <- as.vector(stats::na.omit(y))
#
#   # After removing missing values, both samples need at lest length 5
#   if (length(x) < 5 || length(y) < 5) {
#     stop("Both samples need at least 5 non-missing values.")
#   }
#
#   return(list(x = x, y = y))
# }

#' @title Select principle for computing null distribution
#'
#' @description
#' \code{select_method} is a helper function that chooses the principle for
#' computing the null distribution of a two-sample test.
#'
#' @template x
#' @template y
#' @template test_name
#'
#' @details
#' When the principle is specified by the user, i.e. method contains only
#' one element, the selected method is returned. Otherwisse, if the user
#' does not specify the principle, it depends on the sample size: When both
#' samples contain more than 30 observations, an asymptotic test is perfomed.
#' If one of the samples contains less than 30 observations, the null
#' distribution is computed via the randomization principle. The number of
#' replications \code{n.rep} for the randomization test needs to be specified
#' outside of this function. Each test function contains the argument
#' \code{n.rep} where this can be done.
#'
#' @return A character string which contains the principle for computing the
#'         null distribution.
#'
#' @keywords internal

select_method <- function(x, y, method, test.name) {

  ## Select principle for computing the null distribution ----
  ## If 'method' contains only one entry, the principle is specified by the
  ## user
  if (length(method) == 1) {
    # User-specified principle
    return(method)
  }

  if (test.name %in% c("hl1_test", "hl2_test", "med_test", "trimmed_test", "m_test")) {
    # Automatic selection of the principle, if not specified by the user
    if (length(method) > 1 & identical(method, c("asymptotic", "permutation", "randomization"))) {
      if (length(x) >= 30 & length(y) >= 30) {
        method <- "asymptotic"
      } else {
        method <- "randomization"
      }
    }
    return(method)
  } else {
    stop(paste("Automatic selection not implemented for test",
               paste0("'", test.name, "'."),
               "Please specify 'method' explictly."))
  }
}

#' @title Finite-sample test decision for HL1-, HL2-, and MED-tests
#'
#' @description
#' \code{compute_results_finite} is a helper function to compute the test
#' decision for the HL1-, HL2-, and MED-test, when \code{method = "randomization"}
#' or \code{method = "permutation"}.
#'
#' @template x
#' @template y
#' @template alternative
#' @template delta
#' @template method
#' @template n_rep
#' @template type_rob_perm
#'
#' @return A list containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{estimates}{the location estimates for both samples in case of the HL1-
#'                 and the MED-tests. The estimate for the location difference
#'                 in case of the HL2-tests.}
#' \item{p.value}{the p-value for the test.}
#'
#' @keywords interal

compute_results_finite <- function(x, y, alternative, delta, method, n.rep, type) {

  ## Compute value of the test statistic ----
  y <- ifelse(alternative == "two.sided")
  perm.stats <- rob_perm_statistic(x, y + delta, type = type, na.rm = TRUE)
  statistic <- perm.stats$statistic

  ## Compute location estimates ----
  # For HL1- and MED-tests: location estimates
  # For HL2-tests: Estimator for location difference
  if (type %in% c("HL11", "HL12", "MED1", "MED2")) {
    estimates <- perm.stats$estimates
    estimates[2] <- estimates[2] - delta
  } else if (type %in% c("HL21", "HL22")) {
    estimates <- hodges_lehmann_2sample(x, y)
  }

  # Compute permutation or randomization distribution ----
  distribution <- suppressWarnings(
    perm_distribution(x = x,
                      y = y + delta,
                      type = type,
                      randomization = (method == "randomization"),
                      n.rep = n.rep
    )
  )

  ## Compute p-value ----
  p.value <- calc_perm_p_value(
    statistic,
    distribution,
    m = length(x),
    n = length(y),
    randomization = (method == "randomization"),
    n.rep = n.rep,
    alternative = alternative
  )

  return(list(statistic = statistic, estimates = estimates, p.value = p.value))
}

#' @title Test decision for asymptotic versions of HL1-, HL2-, and MED-tests
#'
#' @description
#' \code{compute_results_asymptotic} is a helper function to compute the test
#' decision for the HL1-, HL2-, and MED-test, when \code{method = "asymptotic"}.
#'
#' @template x
#' @template y
#' @template alternative
#' @template delta
#' @template type_rob_perm
#'
#' @return A list containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{estimates}{the location estimates for both samples in case of the HL1-
#'                 and the MED-tests. The estimate for the location difference
#'                 in case of the HL2-tests.}
#' \item{p.value}{the p-value for the test.}
#'
#' @keywords interal

compute_results_asymptotic <- function(x, y, alternative, delta, type) {

  ## Sample sizes ----
  m <- length(x)
  n <- length(y)
  lambda <- m/(m + n)

  ## Test statistic and estimates for HL1- and HL2-tests ----
  if (type %in% c("HL11", "HL12", "HL21", "HL22")) {
    # Kernel-density estimation for density of pairwise differences
    xcomb <- utils::combn(x, 2)
    ycomb <- utils::combn(y + delta, 2)
    pwdiffs <- c(xcomb[2, ] - xcomb[1, ], ycomb[2, ] - ycomb[1, ])
    dens <- stats::density(pwdiffs)
    dens <- stats::approxfun(dens)

    # Density estimate at 0
    int <- dens(0)

    ## Test statistic and location estimates (HL1)/estimate for location
    ## difference (HL2)
    if (type %in% c("HL11", "HL12")) {
      estimates <- c(hodges_lehmann(x), hodges_lehmann(y + delta))
      statistic <- sqrt(12 * lambda * (1 - lambda)) * int * sqrt(m + n) * (estimates[1] - estimates[2])
      estimates[2] <- estimates[2] - delta
    } else if (type %in% c("HL21", "HL22")) {
      estimates <- hodges_lehmann_2sample(x, y + delta)
      statistic <- sqrt(12 * lambda * (1 - lambda)) * int * sqrt(m + n) * estimates
      estimates <- hodges_lehmann_2sample(x, y)
    }
  } else if (type %in% c("MED1", "MED2")) {
    ## Test statistic and estimates for MED-tests ----

    # Sample medians
    med.x <- stats::median(x, na.rm = TRUE)
    med.y <- stats::median(y + delta, na.rm = TRUE)

    # Kernel-density estimation for differences from sample median
    diff <- c(x - med.x, y + delta - med.y)
    dens <- stats::approxfun(stats::density(diff))
    med <- dens(0)

    # Location estimates and test statistic
    estimates <- c(med.x, med.y - delta)
    statistic <- sqrt(m*n/(m+n)) * 2 * med * (estimates[1] - estimates[2])
  }

  ## Compute p-value ----
  p.value <- switch (alternative,
                     two.sided = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
                     greater = stats::pnorm(statistic, lower.tail = FALSE),
                     less = stats::pnorm(statistic, lower.tail = TRUE)
  )

  return(list(statistic = statistic, estimates = estimates, p.value = p.value))
}
