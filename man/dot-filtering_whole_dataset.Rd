% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtering.R
\name{.filtering_whole_dataset}
\alias{.filtering_whole_dataset}
\title{Filtering data set by presence of values across whole data set}
\usage{
.filtering_whole_dataset(data, threshold)
}
\arguments{
\item{data}{The data frame to run filtering on.}

\item{threshold}{The threshold how many samples need to be present in the
subset the different methods look at to retain a data point.}
}
\value{
A filtered data frame.
}
\description{
Filtering data set by presence of values across whole data set
}
