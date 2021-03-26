#' Preprocess data for the robust two sample tests
#'
#' The function performs all the necessary preprocessing steps for the data.
#'
#' @template x
#' @temlate y
#' @template delta
#' @template na_rm
#' @template wobble
#' @template seed
#' @template var_test
#'
#' @returns A named list containing (the possibly transformed) \code{x}, \code{y} and \code{delta}
#'
#' @keywords internal

preprocess_data <- function(x, y, delta, na.rm, wobble, wobble.seed, var.test) {

  if (na.rm) {
    if (any(is.na(x)) | any(is.na(y))) {
      x <- na.omit(x)
      y <- na.omit(y)
    }
  } else if (!na.rm) {
    if (any(is.na(x)) | any(is.na(y))) {
      return(NA)
    }
  }

  if (wobble) {

    if (!(length(unique(x)) == length(x) &
          length(unique(y)) == length(y) &
          length(unique(c(x, y))) == length(c(x, y)))) {

      if (missing(wobble.seed) | is.null(wobble.seed)) {
        wobble.seed <- sample(1e6, 1)
      }
      set.seed(wobble.seed)


      xy <- wobble(x, y, check = FALSE)
      x <- xy$x
      y <- xy$y

      warning(paste0("Added random noise to x and y. The seed is ",
                     wobble.seed, "."))
    }
  }

  if (var.test) {

    if (any(c(x, y) == 0)) {

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
    ## The transform of Delta needs to be done in the function!!
  }


  return(list(x = x, y = y, delta = delta))
}
## Rückgabe von NA auch in den Test-Funktionen selbst auffangen?



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

  checkmate::assert_choice(test.name, choices = c("hl1_test", "hl2_test", "med_test", "trimmed_test", "m_test"))

  ## Checks that are necessary for all tests ----
  checkmate::assert_numeric(x, finite = TRUE, all.missing = FALSE, min.len = 5, null.ok = FALSE)
  checkmate::assert_numeric(y, finite = TRUE, all.missing = FALSE, min.len = 5, null.ok = FALSE)
  checkmate::assert_subset(alternative, choices = c("two.sided", "greater", "less"), empty.ok = FALSE)
  checkmate::assert_subset(method, choices = c("asymptotic", "permutation", "randomization"), empty.ok = FALSE)
  checkmate::assert_subset(scale, choices = c("S1", "S2"), empty.ok = FALSE)
  checkmate::assert_count(n.rep, na.ok = FALSE, positive = TRUE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)
  checkmate::assert_flag(var.test, na.ok = FALSE, null.ok = FALSE)
  if (var.test) {
    checkmate::assert_numeric(delta, lower = 0, finite = TRUE, any.missing = FALSE, len = 1, null.ok = FALSE)
  } else if (!var.test) {
    checkmate::assert_numeric(delta, finite = TRUE, any.missing = FALSE, len = 1, null.ok = FALSE)
  }
  checkmate::assert_flag(wobble, na.ok = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(wobble.seed, finite = TRUE, any.missing = FALSE, len = 1, null.ok = TRUE)

  ## Additional checks for trimmed t-test and M-tests ---

  # Trimmed t-test
  if (test.name == "trimmed_test") {
    checkmate::assert_number(gamma, na.ok = FALSE, lower = 0, upper = 0.5, finite = TRUE, null.ok = FALSE)
  }

  # M-tests
  if (test.name == "m_test") {
    checkmate::assert_subset(psi, choices = c("huber", "hampel", "bisquare"), empty.ok = FALSE)
    checkmate::assert_numeric(k, lower = 0, len = ifelse(psi == "hampel", 3, 1), finite = TRUE, any.missing = FALSE, null.ok = FALSE)
  }
}


## Remove missing values in 'x' and 'y' ----
remove_missing_values <- function(x, y, na.rm) {

  if (!na.rm & (any(is.na(x)) || any(is.na(y)))) {
    return(NA_real_)
  } else if (na.rm & (any(is.na(x)) || any(is.na(y)))) {
    x <- as.vector(stats::na.omit(x))
    y <- as.vector(stats::na.omit(y))

    # After removing missing values, both samples need at lest length 5
    if (length(x) < 5 || length(y) < 5) {
      stop("Both samples need at least 5 non-missing values.")
    }
  }

  return(list(x = x, y = y))

}

select_method <- function(x, y, method, test.name) {
  if (length(method) == 1) {
    return(method)
  }

  if (test.name %in% c("hl1_test", "hl2_test", "med_test", "trimmed_test", "m_test")) {
    if (length(method) > 1 & identical(method, c("asymptotic", "permutation", "randomization"))) {
      if (length(x) >= 30 & length(y) >= 30) {
        method <- "asymptotic"
      } else {
        method <- "randomization"
      }
    }
  }

  return(method)
}

compute_results_finite <- function(x, y, alternative, delta, method, n.rep, type) {

  # Compute values of the test statistic and location estimates for both
  # samples
  perm.stats <- rob_perm_statistic(x, y + delta, type = type, na.rm = TRUE)

  statistic <- perm.stats$statistic

  if (type %in% c("HL11", "HL12", "MED1", "MED2")) {
    estimates <- perm.stats$estimates
    estimates[2] <- estimates[2] - delta
  } else if (type %in% c("HL21", "HL22")) {
    estimates <- hodges_lehmann_2sample(x, y)
  }

  # Compute permutation or randomization distribution
  distribution <- suppressWarnings(
    perm_distribution(x = x,
                      y = y + delta,
                      type = type,
                      randomization = (method == "randomization"),
                      n.rep = n.rep
    )
  )

  ## Compute p-value
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

compute_results_asymptotic <- function(x, y, alternative, delta, type) {
  # Sample sizes
  m <- length(x)
  n <- length(y)
  lambda <- m/(m + n)

  if (type %in% c("HL11", "HL12", "HL21", "HL22")) {
    # Kernel-density estimation for density of pairwise differences
    xcomb <- utils::combn(x, 2)
    ycomb <- utils::combn(y + delta, 2)
    pwdiffs <- c(xcomb[2, ] - xcomb[1, ], ycomb[2, ] - ycomb[1, ])
    dens <- stats::density(pwdiffs)
    dens <- stats::approxfun(dens)

    int <- dens(0)

    # Compute values of the test statistic and location estimates for both
    # samples
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
    med.x <- stats::median(x, na.rm = TRUE)
    med.y <- stats::median(y + delta, na.rm = TRUE)

    diff <- c(x - med.x, y + delta - med.y)
    dens <- stats::approxfun(stats::density(diff))
    med <- dens(0)

    estimates <- c(med.x, med.y - delta)

    statistic <- sqrt(m*n/(m+n)) * 2 * med * (estimates[1] - estimates[2])
  }

  ## Compute p-value
  p.value <- switch (alternative,
                     two.sided = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
                     greater = stats::pnorm(statistic, lower.tail = FALSE),
                     less = stats::pnorm(statistic, lower.tail = TRUE)
  )

  return(list(statistic = statistic, estimates = estimates, p.value = p.value))
}

