% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testStatistics.R
\name{rob_perm_statistic}
\alias{rob_perm_statistic}
\title{Robust test statistics based on robust location estimators}
\usage{
rob_perm_statistic(
  x,
  y,
  type = c("HL11", "HL12", "HL21", "HL22", "MED1", "MED2"),
  na.rm = FALSE
)
}
\arguments{
\item{x}{a (non-empty) numeric vector of data values.}

\item{y}{a (non-empty) numeric vector of data values.}

\item{type}{a character string specifying the desired test statistic. It must
be one of \code{"HL11"} (default),
\code{"HL12", "HL21", "HL22", "MED1"}, and \code{"MED2"}, where
\code{"HL1"}, \code{"HL2"} and \code{"MED"} specify the location
estimator and the numbers \code{1} and \code{2} the scale
estimator, see the vignette
\code{vignette("robnptests")} for more information.}

\item{na.rm}{a logical value indicating whether NA values in \code{x} and
\code{y} should be stripped before the computation proceeds. The
default is \code{na.rm = FALSE}.}
}
\value{
A named list containing the following components:
        \item{statistic}{the selected test statistic.}
        \item{estimates}{estimate of location for each sample if available.}
}
\description{
\code{rob_perm_statistic} calculates test statistics for robust
permutation/randomization tests based on the sample median, the one-sample
Hodges-Lehmann estimator, or the two-sample Hodges-Lehmann estimator.
}
\details{
The test statistics returned by \code{rob_perm_statistic} are of the
         form \deqn{D_i/S_j} where the D_i, i = 1,...,3, are different
         estimators of location and the S_j, j = 1,...,4, are estimates for
         the mutual sample scale. See \insertCite{FriDeh11robu;textual}{robnptests}
         or the vignette \code{vignette("robnptests")} for details.
}
\examples{
# Generate random samples
set.seed(108)
x <- rnorm(20); y <- rnorm(20)

# Compute HL21-statistic
rob_perm_statistic(x, y, type = "HL21")

}
\references{
\insertRef{FriDeh11robu}{robnptests}
}
