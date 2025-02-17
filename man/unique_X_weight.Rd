% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iclogcondist_functions.R
\name{unique_X_weight}
\alias{unique_X_weight}
\title{Find Unique Rows in a Matrix and Their Weights}
\usage{
unique_X_weight(X)
}
\arguments{
\item{X}{A matrix. The matrix whose unique rows are to be found.}
}
\value{
A list containing two components:
\describe{
  \item{unique_X}{A matrix of the unique rows from the input matrix.}
  \item{weight}{An integer vector containing the frequency (weight) of each unique row.}
}
}
\description{
This function finds the unique rows of a given matrix and calculates the frequency (weight) of each unique row. 
It returns both the unique rows and the weights (the number of occurrences of each row).
}
