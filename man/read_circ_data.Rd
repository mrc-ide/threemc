% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{read_circ_data}
\alias{read_circ_data}
\title{Function to read in Circumcision Data}
\usage{
read_circ_data(path, filters = NULL, selected = NULL, ...)
}
\arguments{
\item{path}{Path to data.}

\item{filters}{Optional named vector, whose values dictate the values
filtered for in the corresponding column names. Only supports filtering for
one value for each column. default: NULL}

\item{selected}{Optional columns to select, removing others, default = NULL}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
relevant data set, filtered as desired.
}
\description{
Function to read in circumcision data to fit model. Handles
csv with \code{\link[data.table]{fread}} (but outputs data as a
\code{data.frame}), and geographical data with code{\link[sf]{read_sf}} (for which
it also adds unique identifiers for each \code{area_level}).
}
\seealso{
\code{\link[data.table]{fread}}
\code{\link[sf]{read_sf}}
}
