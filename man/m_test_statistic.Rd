% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testStatistics.R
\name{m_test_statistic}
\alias{m_test_statistic}
\title{Test statistics for the M-tests}
\usage{
m_test_statistic(x, y, psi, k = robustbase::.Mpsi.tuning.default(psi), ...)
}
\arguments{
\item{x}{a (non-empty) numeric vector of data values.}

\item{y}{a (non-empty) numeric vector of data values.}

\item{psi}{kernel used for optimization.
Must be one of \code{"bisquare"}, \code{"hampel"} and \code{"huber"}.
The default is \code{"huber"}.}

\item{k}{tuning parameter(s) for the respective kernel function,
defaults to parameters implemented in \code{\link[robustbase:lmrob.control]{.Mpsi.tuning.default(psi)}}
in the package \href{https://cran.r-project.org/package=robustbase}{robustbase}.}

\item{...}{additional arguments \code{c1} and \code{c2} that can be passed to
the function \code{scaleTau2()}, which is used internally for
estimating the within-sample dispersion, in order to account for
non-normal distributions; see \insertCite{MarZam02robu;textual}{robnptests}.}
}
\value{
A named list containing the following components:
        \item{statistic}{standardized test statistic.}
        \item{estimates}{M-estimates of location for both \code{x} and \code{y}.}
}
\description{
\code{m_test_statistic} calculates the test statistics for
tests based on M-estimators.
}
\details{
For details on how the test statistic is constructed, we refer to the
vignette \code{vignette("m_tests")}
}
\examples{
# Generate random samples
set.seed(108)
x <- rnorm(20); y <- rnorm(20)

# Compute Huber-M-statistic
m_test_statistic(x, y, psi = "huber")

}
