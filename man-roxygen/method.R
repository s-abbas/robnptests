#' @param method a character string specifying how the p-value is computed,
#'               \code{"asymptotic"} for an asymptotic test based on a normal
#'               approximation, \code{"permutation"} for a permutation test and
#'               \code{"randomization"} for a randomization test.
#'               The permutation test uses all splits of the joint sample into two samples
#'               of sizes \code{m} and \code{n}, while the randomization test draws
#'               \code{n.rep} random splits with replacement.
