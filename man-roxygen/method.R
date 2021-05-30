#' @param method a character string specifying how the p-value is computed with
#'               possible values  \code{"asymptotic"} for an asymptotic test
#'               based on a normal approximation, \code{"permutation"} for a
#'               permutation test, and \code{"randomization"} for a randomization
#'               test. The permutation test uses all splits of the joint sample
#'               into two samples of sizes \code{m} and \code{n}, while the
#'               randomization test draws \code{n.rep} random splits with
#'               replacement. The values \code{m} and \code{n} denote the
#'               sample sizes.
#'               If not specified explicitly, defaults to
#'               \code{"permutation"}  if \code{m < 30}, \code{n < 30} and
#'                \code{n.rep >= choose(m + n, m)},
#'               \code{"randomization"} if \code{m < 30}, \code{n < 30} and
#'                \code{n.rep < choose(m + n, m)}, and
#'                \code{"asymptotic"} if \code{m >= 30} and \code{n >= 30}.
