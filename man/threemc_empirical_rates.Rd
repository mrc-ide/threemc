% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shell_dataset.R
\name{threemc_empirical_rates}
\alias{threemc_empirical_rates}
\title{Use model shell dataset to estimate empirical circumcision rates}
\usage{
threemc_empirical_rates(
  out,
  areas,
  area_lev,
  populations,
  age_groups = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
    "40-44", "45-49", "50-54", "54-59", "0+", "10+", "15+", "15-24", "10-24", "15-29",
    "10-29", "15-39", "10-39", "15-49", "10-49")
)
}
\arguments{
\item{out}{Shell dataset outputted by \link[threemc]{create_shell_dataset}}

\item{areas}{\code{sf} shapefiles for specific country/region.}

\item{area_lev}{\itemize{
\item Desired admin boundary level to perform the analysis on.
}}

\item{populations}{\code{data.frame} containing populations for each
region in tmb fits.}

\item{age_groups}{Age groups to aggregate by, Default:
c("0-4",   "5-9",   "10-14", "15-19",
"20-24", "25-29", "30-34", "35-39",
"40-44", "45-49", "50-54", "54-59",
"15-24", "10-24", "15-29", "10-29",
"15-39", "10-39", "15-49", "10-49"
)}
}
\description{
Takes the shell dataset with a row for every unique area ID,
area name, year and circumcision age in survey data outputed by
\link[threemc]{create_shell_dataset} and returns the empirical circumcision
rates for each row, aggregated to age groups from single ages. Also converts
from wide format to long format.
}
\seealso{
\code{\link[threemc]{create_shell_dataset}}
}
