#' @param method a character string specifying the test method used: \code{"permutation"}
#'               for a permutation and \code{"randomization"} for a randomization test.
#'               The permutation test uses all data splits into two samples,
#'               while the randomization test draws \code{n.rep} random splits
#'               with replacement. If not specified explicitly,
#'               defaults to
#'               \code{"permutation"},  if \code{m < 30}, \code{n < 30} and
#'                \code{n.rep >= choose(m+n, m)},
#'               \code{"randomization"}, if \code{m < 30}, \code{n < 30} and
#'                \code{n.rep < choose(m+n, m)}, and
#'                \code{"asymptotic"}, if \code{m >= 30} and \code{n >= 30}.

