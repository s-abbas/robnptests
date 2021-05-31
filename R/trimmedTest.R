## ----------------------------------------------------------------------------
## Trimmed t-Test (Yuen's t-test)
## ----------------------------------------------------------------------------

#' @title Two-sample trimmed t-test (Yuen's t-Test)
#'
#' @description
#' \code{trimmed_test} performs the two-sample trimmed t-test.
#'
#' @template x
#' @template y
#' @template gamma_trimmed
#' @template alternative
#' @template delta
#' @template method
#' @template n_rep
#' @template na_rm
#' @template var_test
#' @template wobble_seed
#'
#' @details
#' The function performs Yuen's t-test based on the trimmed mean and winsorized
#' variance \insertCite{YueDix73appr}{robnptests}.
#' The amount of trimming / winsorization is set in \code{gamma} and
#' defaults to 0.2, i.e. 20\% of the values are removed/replaced.
#' In addition to the asymptotic distribution we provide a permutation and a randomization
#' version of the test.
#'
#' When computing a randomization distribution based on randomly drawn splits
#' with replacement, the function \code{\link[statmod]{permp}} \insertCite{PhiSmy10perm}{robnptests}
#' is used to calculate the p-value.
#'
#' For \code{var.test = TRUE}, the test compares the two samples for a difference in scale.
#' This is achieved by log-transforming the original observations so that a potential
#' scale difference appears as a location difference between the transformed samples;
#' see \insertCite{Fri12onli;textual}{robnptests}. The sample should not contain zeros
#' to prevent problems with the necessary log-transformation. If it contains zeros,
#' uniform noise is added to all variables in order to remove zeros. A warning is
#' printed.
#'
#' If the sample has been modified because of zeros when \code{var.test = TRUE},
#' the modified samples can be retrieved using
#'
#' \code{set.seed(wobble.seed); wobble(x, y)}
#'
#' Both samples need to contain at least 5 non-missing values.
#'
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom for the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the trimmed means of \code{x} and \code{y}.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating how the p-value was computed.}
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
#' \insertRef{YueDix73appr}{robnptests}
#'
#' \insertRef{Yue74trim}{robnptests}
#'
#' \insertRef{Fri12onli}{robnptests}
#'
#' @export

trimmed_test <- function(x, y, gamma = 0.2,
                         alternative = c("two.sided", "less", "greater"),
                         method = c("asymptotic", "permutation", "randomization"),
                         delta = ifelse(var.test, 1, 0),
                         n.rep = 1000,
                         na.rm = FALSE, var.test = FALSE,
                         wobble.seed = NULL) {

  # Check input arguments ----
  check_test_input(x = x, y = y, gamma = gamma, alternative = alternative,
                   method = method, delta = delta, n.rep = n.rep, na.rm = na.rm,
                   var.test = var.test, wobble = FALSE, wobble.seed = wobble.seed,
                   test.name = "trimmed_test")

  # Extract names of data sets ----
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  # Match 'alternative' ----
  # 'method' not matched because computation of p-value depends on sample sizes
  # if no value is specified by the user
  alternative <- match.arg(alternative)

  # Data preprocessing ----
  prep <- preprocess_data(x = x, y = y, delta = delta, na.rm = na.rm,
                          wobble = FALSE, wobble.seed = wobble.seed,
                          var.test = var.test)

  if (!all(is.na(prep))) {
    x <- prep$x
    y <- prep$y
    delta <- prep$delta
  } else {
    return(NA)
  }

  # Select method for computing the p-value ----
  method <- select_method(x = x, y = y, method = method, test.name = "trimmed_test",
                          n.rep = n.rep)

  # Test decision ----

  # Test statistic, location estimates for both samples, and degrees of freedom
  t.stats <- trimmed_t(x, y + delta, gamma = gamma, na.rm = na.rm)
  statistic <- t.stats$statistic
  estimates <- t.stats$estimates
  estimates[2] <- t.stats$estimates[2] - delta
  df <- t.stats$df

  if (method %in% c("permutation", "randomization")) {
    # Test decision for permutation or randomization test

    m <- length(x)
    n <- length(y)
    y <- y + delta
    complete <- c(x, y)

    if (method == "permutation") {
      # Compute permutation distribution
      splits <- gtools::combinations((m + n), m, 1:(m + n))

      distribution <- apply(splits, 1, function(s) {
        trimmed_t(x = complete[s], y = complete[-s], gamma = gamma)$statistic
      })
    } else if (method == "randomization") {
      # Compute randomization distribution
      n.rep <- min(choose(length(x) + length(y), length(x)), n.rep)
      splits <- replicate(n.rep, sample(complete))

      distribution <- apply(splits, 2, function(s) {
        trimmed_t(x = s[1:m], y = s[(m + 1):(m + n)], gamma = gamma)$statistic
      })
    }

    p.value <- calc_perm_p_value(statistic, distribution, m = m, n = n,
                                 randomization = (method == "randomization"),
                                 n.rep = n.rep, alternative = alternative)

  } else if (method == "asymptotic") {
    # Test decision for asymptotic test
    p.value <- switch(alternative,
                      two.sided = 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE),
                      greater = stats::pt(statistic, df = df, lower.tail = FALSE),
                      less = stats::pt(statistic, df = df, lower.tail = TRUE)
    )
  }

  # Prepare output ----

  # Assign names to results
  if (var.test) {
    names(estimates) <- c("Trimmed mean of log(x^2)", "Trimmed mean of log(y^2)")
    names(delta) <- "ratio of squared scale parameters"
    delta <- exp(delta)
  } else {
    names(estimates) <- c("Trimmed mean of x", "Trimmed mean of y")
    names(delta) <- "location shift"
  }

  names(statistic) <- "trimmed t"
  names(df) <- "df"

  # Information on applied test
  if (method == "randomization") {
    method <- paste0("Randomization test based on trimmed means ", "(", n.rep, " random permutations)")
  } else if (method == "permutation") {
    method <- "Exact permutation test based on trimmed means"
  } else {
    method <- "Yuen's t-test"
  }

  # Results
  res <- list(statistic = statistic, parameter = df, p.value = p.value,
              estimate = estimates, null.value = delta, alternative = alternative,
              method = method, data.name = dname)
  class(res) <- "htest"

  return(res)
}
