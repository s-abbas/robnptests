#' @param wobble logical value indicating whether uniform noise should be added to the sample before
#'               proceeding with the test, in order to make sure there are no duplicate values in the
#'               sample that can cause the scale estimate to be 0. Only necessary for the permutation and
#'               randomization version of the test. The default is \code{wobble = FALSE}.
#' @param wobble.seed integer value used as a seed for the random number generation in case of \code{wobble = TRUE}.
#'                    If no seed is given, it is set randomly. In any case, the
#'                    seed used will be returned in a warning message. The argument is
#'                    ignored if wobble = FALSE.
