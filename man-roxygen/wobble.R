#' @param wobble logical value indicating whether the sample should be checked for
#'               duplicated values that can cause the scale estimate to be zero.
#'               If such values are present, uniform noise is added to the sample,
#'               see \code{\link[robTests]{wobble()}}.
#'               Only necessary for the permutation and randomization version of the test.
#'               The default is \code{wobble = FALSE}.
