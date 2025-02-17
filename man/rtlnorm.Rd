% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iclogcondist_simulation_dist.R
\name{rtlnorm}
\alias{rtlnorm}
\title{Simulate from a Truncated Log-Normal Distribution}
\usage{
rtlnorm(n, meanlog = 0, sdlog = 1, upper_bound = Inf)
}
\arguments{
\item{n}{An integer specifying the number of random samples to generate.}

\item{meanlog}{A numeric value representing the mean of the log-normal distribution on the log scale. Default is \code{0}.}

\item{sdlog}{A positive numeric value representing the standard deviation of the log-normal distribution on the log scale. Default is \code{1}.}

\item{upper_bound}{A positive numeric value indicating the upper truncation point. Default is \code{Inf} (no truncation).}
}
\value{
A numeric vector of \code{n} random samples from the truncated log-normal distribution.
}
\description{
This function generates random samples from a truncated log-normal distribution
using an acceptance-rejection method.
}
\examples{
# Generate 10 random samples from a truncated log-normal distribution
rtlnorm(10, meanlog = 0, sdlog = 1, upper_bound = 5)

}
