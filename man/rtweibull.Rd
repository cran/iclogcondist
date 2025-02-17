% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iclogcondist_simulation_dist.R
\name{rtweibull}
\alias{rtweibull}
\title{Simulate from a Truncated Weibull Distribution}
\usage{
rtweibull(n, shape = 1, scale = 1, upper_bound = Inf)
}
\arguments{
\item{n}{An integer specifying the number of random samples to generate.}

\item{shape}{A positive numeric value representing the shape parameter of the Weibull distribution. Default is \code{1}.}

\item{scale}{A positive numeric value representing the scale parameter of the Weibull distribution. Default is \code{1}.}

\item{upper_bound}{A positive numeric value indicating the upper truncation point. Default is \code{Inf} (no truncation).}
}
\value{
A numeric vector of \code{n} random samples from the truncated Weibull distribution.
}
\description{
This function generates random samples from a truncated Weibull distribution
using inverse transform sampling. When \code{shape = 1}, it reduces to a truncated exponential distribution.
}
\examples{
# Generate 10 random samples from a truncated Weibull distribution
rtweibull(10, shape = 2, scale = 1, upper_bound = 5)

}
