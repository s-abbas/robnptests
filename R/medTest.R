## ----------------------------------------------------------------------------
## Two sample test based on median differences
## ----------------------------------------------------------------------------

#' Two-sample location tests based on the sample median
#'
#' @description
#' \code{med_test} performs a two-sample location test based on
#' the difference of the sample medians for both samples.
#'
#' @template x
#' @template y
#' @template alternative
#' @template delta
#' @template method
#' @template scale_med
#' @template n_rep
#' @template na_rm
#' @template var_test
#' @template wobble
#' @template seed
#'
#' @details
#' When computing a randomization distribution based on randomly drawn splits with replacement, the results of
#' Smyth & Phipson (2010) to calculate the p-value are used. The test statistics and the asymptotic distribution are taken from Fried & Dehling (2011).
#'
#' The estimator of the location shift is the difference of the medians of \code{x} and \code{y}.
#'
#' The test statistics for the permutation and randomization version of the test is standardized using a robust scale estimator.
#' \code{scale = "S3"} represents use of
#'
#' \deqn{S = 2 * ( |X_1 - med(X)|,...,|X_m - med(X)|, |Y_1 - med(Y)|,...,|Y_n - med(Y)| ),}
#'
#' \code{scale = "S4"} uses
#'
#' \deqn{S = ( med( |X_1 - med(X)|,...,|X_m - med(X)| ) + med( |Y_1 - med(Y)|,...,|Y_n - med(Y)| ). }
#'
#' For more details see Fried & Dehling (2011).
#'
#' For \code{var.test = TRUE}, the test compares the two samples for a difference in scale.
#' This is achieved by log-transforming the original observations so that a potential
#' scale difference appears as a location difference between the transformed samples;
#' see Fried (2012). The sample cannot contain zeros due to the necessary log-transformation.
#' If it contains zeros, uniform noise is added to all variables in order to remove zeros.
#' A warning is printed.
#'
#' If the sample has been modified (either because of 0's for \code{var.test = TRUE}, or
#' as \code{wobble = TRUE}, the modified samples can be retrieved using
#'
#' \code{set.seed(wobble.seed); wobble(x, y)}
#'
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the sample medians of \code{x} and \code{y}.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of test was performed.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{PhiSmy10perm}{robTests}
#'
#' \insertRef{FriDeh11robu}{robTests}
#'
#' @examples
#' ## Generate random samples
#' set.seed(108)
#' x <- rnorm(20); y <- rnorm(20)
#'
#' ## Asymptotic MED test
#' med_test(x, y, method = "asymptotic", scale = "S3")
#'
#' ## MED2 test using randomization principle by drawing 1000 random permutations
#' ## with replacement
#'
#' \dontrun{
#' med_test(x, y, method = "randomization", n.rep = 1000, scale = "S4")
#' }
#'
#' @export

med_test <- function(x, y, alternative = c("two.sided", "greater", "less"),
                     delta = ifelse(var.test, 1, 0),
                     method = c("asymptotic", "permutation", "randomization"),
                     scale = c("S3", "S4"), n.rep = 10000,
                     na.rm = FALSE, var.test = FALSE,
                     wobble = FALSE, wobble.seed = NULL) {

  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  if (!na.rm & (any(is.na(x)) | any(is.na(y)))) {
    return(NA)
  } else if (na.rm & (any(is.na(x)) | any(is.na(y)))) {
    x <- as.numeric(stats::na.omit(x))
    y <- as.numeric(stats::na.omit(y))
  }

  if (length(x) < 5 || length(y) < 5) {
    stop("Both samples need at least 5 non-missing values.")
  }


  if (wobble) {

    if (is.null(wobble.seed)) wobble.seed <- sample(1e6, 1)
    set.seed(wobble.seed)

    xy <- wobble(x, y)
    x <- xy$x
    y <- xy$y

    warning(paste0("Added random noise to x and y. The seed is ",
                   wobble.seed, "."))
  }

  ## If necessary: Transformation to test for difference in scale
  if (var.test) {

    if (any(c(x, y) == 0)) {

      if (is.null(wobble.seed)) wobble.seed <- sample(1e6, 1)
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
  }

  alternative <- match.arg(alternative)
  # method <- match.arg(method)
  scale <- match.arg(scale)

  if (scale == "S3") {
    type <- "MED1"
  } else if (scale == "S4") {
    type <- "MED2"
    } else stop(" 'scale' must one of 'S3' and 'S4' ")

  ## Error handling
  if (!missing(delta) && (length(delta) != 1 || is.na(delta))) {
    stop ("'delta' must be a single number.")
  }

  if (length(method) > 1 & identical(method, c("asymptotic", "permutation", "randomization"))) {
    if (length(x) >= 30 & length(y) >= 30) method <- "asymptotic"
    else method <- "randomization"
  }

  if (!(method %in% c("asymptotic", "permutation", "randomization"))) {
    stop (" 'method' must be one of 'asymptotic', 'permutation' or 'randomization'. ")
  }

  if (method %in% c("permutation", "randomization")) {

    ## Results of rob_perm_statistic
    perm.stats <- rob_perm_statistic(x, y + delta, type = type, na.rm = na.rm)

    statistic <- perm.stats$statistic
    estimates <- perm.stats$estimates

    if (delta != 0) estimates[2] <- stats::median(y)

    ## Calculate permutation distribution
    # if (method == "randomization") {
    #   randomization <- TRUE
    # } else {
    #   randomization <- FALSE
    # }

    distribution <- suppressWarnings(perm_distribution(x = x, y = y + delta, type = type,
                                      randomization = (method == "randomization"),
                                      n.rep = n.rep))

    ## p-value
    p.value <- calc_perm_p_value(statistic, distribution, m = length(x),
                                 n = length(y), randomization = (method == "randomization"), n.rep = n.rep,
                                 alternative = alternative)

  } else if (method == "asymptotic") {

    med.x <- stats::median(x, na.rm = na.rm)
    med.y <- stats::median(y + delta, na.rm = na.rm)

    diff <- c(x - med.x, y + delta - med.y)

    dens <- stats::approxfun(stats::density(diff))
    med <- dens(0)

    m <- length(x)
    n <- length(y)

    if (delta != 0) estimates <- c(med.x, stats::median(y)) else estimates <- c(med.x, med.y)
    est <- med.x - med.y

    statistic <- sqrt(m*n/(m+n)) * 2 * med * est

    p.value <- switch (alternative,
                       two.sided = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
                       greater = stats::pnorm(statistic, lower.tail = FALSE),
                       less = stats::pnorm(statistic, lower.tail = TRUE)
    )

  }

  ## Assign names to results

  if (var.test) {
    names(estimates) <- c("Median of log(x^2)", "Median of log(y^2)")
    names(delta) <- "ratio of variances"
    delta <- exp(delta)
  } else {
    names(estimates) <- c("Median of x", "Median of y")
    names(delta) <- "location shift"
  }

  names(statistic) <- ifelse(var.test, "S", "D")

  if (method == "randomization") {
    method = "Randomization test based on sample medians"
  } else if (method == "permutation") {
    method = "Exact permutation test based on sample medians"
  } else method = "Asymptotic test based on sample medians"


  res <- list(statistic = statistic, parameter = NULL, p.value = p.value,
              estimate = estimates, null.value = delta, alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}

