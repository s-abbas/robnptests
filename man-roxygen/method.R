#' @param method a character string specifying the test method used,
#'               \code{"asymptotic"} for an asymptotic test based on the
#'               normal distribution, \code{"permutation"} for a permutation and
#'               \code{"randomization"} for a randomization test.
#'               The permutation test uses all data splits into two samples,
#'               while the randomization test draws \code{n.rep} random splits
#'               with replacement.
