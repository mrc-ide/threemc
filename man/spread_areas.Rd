% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{spread_areas}
\alias{spread_areas}
\title{Spread area hierarchy to wide format}
\usage{
spread_areas(
  areas,
  min_level = min(areas$area_level),
  max_level = max(areas$area_level),
  space = TRUE
)
}
\arguments{
\item{areas}{\code{sf} shapefiles for specific country/region.}

\item{min_level}{integer specifying the minimum level wanted}

\item{max_level}{integer specifying the maximum level wanted}

\item{space}{whether to include "space" columns. Excluding these returns the
same object as \code{naomi::spread_areas}, Default: TRUE}
}
\description{
Spread area hierarchy to wide format
}
\keyword{internal}
