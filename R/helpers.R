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

