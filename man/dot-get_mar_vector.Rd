% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputation.R
\name{.get_mar_vector}
\alias{.get_mar_vector}
\title{Obtain vector of MAR rows.}
\usage{
.get_mar_vector(mx, groups)
}
\arguments{
\item{mx}{A expression matrix with missing values.}

\item{groups}{A vector indicating which columns belong to the same group.}
}
\value{
A logical of length nrow(mx). TRUE values are considered MAR.
}
\description{
Missing at random values are defined as observations where not all values
for a specific group/condition are missing, but only some.
}
