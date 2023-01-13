## ----------------------------------------------------------------------------
## Two-sample Hodges-Lehmann test
## ----------------------------------------------------------------------------

#' @title Two-sample location tests based on two-sample Hodges-Lehmann estimator.
#'
#' @description
#' \code{hl2_test} performs a two-sample location test based on the two-sample
#' Hodges-Lehmann estimator for shift.
#'
#' @template x
#' @template y
#' @template alternative
#' @template delta
#' @template method
#' @template scale_hl
#' @template n_rep
#' @template na_rm
#' @template var_test
#' @template wobble
#' @template wobble_seed
#'
#' @details
#' The test statistic for this test is based on the two-sample Hodges-Lehmann
#' estimator of \code{x} and \code{y}, see
#' \code{\link[robnptests]{hodges_lehmann_2sample}}. Three versions of the test
#' are implemented: randomization, permutation, and asymptotic.
#'
#' The test statistic for the permutation and randomization version of the test
#' is standardized using a robust scale estimator, see
#' \insertCite{FriDeh11robu}{robnptests}.
#'
#' With \code{scale = "S1"}, the scale is estimated by
#'
#' \deqn{S = med(|x_i - x_j|: 1 \le i < j \le m, |y_i - y_j|, 1 \le i < j \le n),}
#'
#' whereas \code{scale = "S2"} uses
#'
#' \deqn{S = med(|z_i - z_j|: 1 \le i < j \le m + n).}
#'
#' Here, \eqn{z = (z_1, ..., z_{m + n}) = (x_1 - med(x), ..., x_m - med(x), y_1 - med(y), ..., y_n - med(y))}
#' is the median-corrected sample.
#'
#' The randomization distribution is based on randomly drawn splits with
#' replacement. The function \code{\link[statmod]{permp}} \insertCite{PhiSmy10perm}{robnptests}
#' is used to calculate the p-value. For the asymptotic test, a transformed version
#' of the HL2-estimator, which asymptotically follows a normal distribution, is
#' used. For more details on the asymptotic test, see \insertCite{FriDeh11robu;textual}{robnptests}.
#'
#' For \code{scale_test = TRUE}, the test compares the two samples for a difference
#' in scale. This is achieved by log-transforming the original squared observations,
#' i.e. \code{x} is replaced by \code{log(x^2)} and \code{y} by \code{log(y^2)}.
#' A potential scale difference then appears as a location difference between
#' the transformed samples, see \insertCite{Fri12onli;textual}{robnptests}.
#' Note, that the samples need to be centred around zero. The sample should not
#' contain zeros to prevent problems with the necessary log-transformation. If
#' it contains zeros, uniform noise is added to all variables in order to remove
#' zeros and a message is printed.
#'
#' If the sample has been modified (either because of zeros if \code{scale_test = TRUE}
#' or \code{wobble = TRUE}), the modified samples can be retrieved using
#'
#' \code{set.seed(wobble.seed); wobble(x, y)}.
#'
#' Both samples need to contain at least 5 non-missing values.
#'
#' @return
#' A named list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the estimated location difference between \code{x} and \code{y}
#'                 (if \code{scale_test = FALSE}) or of \code{log(x^2)} and
#'                 \code{log(y^2)} (if \code{scale_test = TRUE}) based on the
#'                 two-sample Hodges-Lehmann estimator.}
#' \item{null.value}{the specified hypothesized value of the mean difference/squared
#'                   scale ratio.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating how the p-value was computed.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{PhiSmy10perm}{robnptests}
#'
#' \insertRef{FriDeh11robu}{robnptests}
#'
#' \insertRef{Fri12onli}{robnptests}
#'
#' @examples
#' # Generate random samples
#' set.seed(108)
#' x <- rnorm(20); y <- rnorm(20)
#'
#' # Asymptotic HL2 test
#' hl2_test(x, y, method = "asymptotic", scale = "S1")
#'
#' \dontrun{
#' # HL22 test using randomization principle by drawing 1000 random permutations
#' # with replacement
#'
#' hl2_test(x, y, method = "randomization", n.rep = 1000, scale = "S2")
#' }
#'
#' @export

hl2_test <- function(x,
                     y,
                     alternative = c("two.sided", "greater", "less"),
                     delta = ifelse(scale_test, 1, 0),
                     method = c("asymptotic", "permutation", "randomization"),
                     scale = c("S1", "S2"),
                     n.rep = 10000,
                     na.rm = FALSE,
                     scale_test = FALSE,
                     wobble = FALSE,
                     wobble.seed = NULL) {

  # Check input arguments ----
  check_test_input(x = x, y = y, alternative = alternative, delta = delta,
                   method = method, scale = scale, n.rep = n.rep, na.rm = na.rm,
                   scale_test = scale_test, wobble = wobble, wobble.seed = wobble.seed,
                   test.name = "hl2_test")

  # Extract names of data sets ----
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  # Match 'alternative' and 'scale' ----
  # 'method' not matched because computation of p-value depends on sample sizes
  # if no value is specified by the user
  alternative <- match.arg(alternative)
  scale <- match.arg(scale)

  # Data preprocessing ----
  prep <- preprocess_data(x = x, y = y, delta = delta, na.rm = na.rm,
                          wobble = wobble, wobble.seed = wobble.seed,
                          scale_test = scale_test)

  if (!all(is.na(prep))) {
    x <- prep$x
    y <- prep$y
    delta <- prep$delta
  } else {
    return(NA)
  }

  # Select scale estimator ----
  if (scale == "S1") {
    type <- "HL21"
  } else if (scale == "S2") {
    type <- "HL22"
  }

  # Select method for computing the p-value ----
  method <- select_method(x = x, y = y, method = method, test.name = "hl2_test",
                          n.rep = n.rep)

  # Test decision ----
  if (method %in% c("permutation", "randomization")) {
    # Test decision for permutation or randomization test
    n.rep <- min(choose(length(x) + length(y), length(x)), n.rep)
    test.results <- compute_results_finite(x = x, y = y, alternative = alternative,
                                           delta = delta, method = method, type = type,
                                           n.rep = n.rep)

  } else if (method == "asymptotic") {
    # Test decision for asymptotic test
    test.results <- compute_results_asymptotic(x = x, y = y, alternative = alternative,
                                               delta = delta, type = type)
  }

  # Test statistic, location estimates for both samples, and p-value
  statistic <- test.results$statistic
  estimates <- test.results$estimates
  p.value   <- test.results$p.value

  # Prepare output ----

  # Assign names to results
  if (scale_test) {
    names(estimates) <- c("HL2 of log(x^2) and log(y^2)")
    names(delta) <- "ratio of squared scale parameters"
    delta <- exp(delta)
  } else {
    names(estimates) <- c("HL2 of x and y")
    names(delta) <- "location shift"
  }

  names(statistic) <- ifelse(scale_test, "S", "D")

  # Information on applied test
  if (method == "randomization") {
    method <- paste0("Randomization test based on HL2-estimator ", "(", n.rep, " random permutations)")
  } else if (method == "permutation") {
    method <- "Exact permutation test based on HL2-estimator"
  } else method <- "Asymptotic test based on HL2-estimator"

  # Results
  res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
              estimate = estimates, null.value = delta, alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}
