% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{preprocess_data}
\alias{preprocess_data}
\title{Preprocess data for the robust two sample tests}
\usage{
preprocess_data(x, y, delta, na.rm, wobble, wobble.seed, scale.test)
}
\arguments{
\item{x}{a (non-empty) numeric vector of data values.}

\item{y}{a (non-empty) numeric vector of data values.}

\item{delta}{a numeric value indicating the true difference in the location
or scale parameter, depending on whether the test should be
performed for a difference in location or in scale.}

\item{na.rm}{a logical value indicating whether NA values in \code{x} and \code{y}
should be stripped before the computation proceeds.}

\item{wobble}{a logical value indicating whether the sample should be checked for
duplicated values that can cause the scale estimate to be zero.
If such values are present, uniform noise is added to the sample,
see \code{\link[robnptests]{wobble}}.}

\item{wobble.seed}{an integer value used as a seed for the random number
generation in case of \code{wobble = TRUE} or when
\code{scale.test = TRUE} with one of the vectors \code{x} and
\code{y} containing zeros. When no seed is specified, it
is chosen randomly and printed in a message. The argument
is ignored if \code{scale.test = FALSE} and/or \code{wobble = FALSE}.}

\item{scale.test}{a logical value to specify if the samples should be compared
for a difference in scale.}
}
\value{
A named list containing the following components:
        \item{\code{x}}{the (possibly transformed) input vector \code{x}.}
        \item{\code{y}}{the (possibly transformed) input vector \code{y}.}
        \item{\code{delta}}{the (possibly transformed) input value
                            \code{delta}.}
}
\description{
\code{preprocess_data} is a helper function that performs several
preprocessing steps on the data before performing the two-sample tests.
}
\details{
The preprocessing steps include the removal of missing values and, if
specified, wobbling and a transformation of the observations to test for
differences in scale.
}
\keyword{internal}
