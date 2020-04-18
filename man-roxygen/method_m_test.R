#' @param method a character string specifying the test method used: \code{"exact"}
#'               for an exact and \code{"sampled"} for a randomized permutation test.
#'               The exact permutation test uses all data splits into two samples,
#'               while the randomized test draws \code{n.rep} random splits
#'               with replacement.
