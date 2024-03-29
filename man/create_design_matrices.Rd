% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_prep.R
\name{create_design_matrices}
\alias{create_design_matrices}
\title{Create Design Matrices}
\usage{
create_design_matrices(
  dat,
  area_lev = NULL,
  k_dt_age = 5,
  k_dt_time = 5,
  inc_time_tmc = FALSE
)
}
\arguments{
\item{dat}{Shell dataset (datputted by \link[threemc]{create_shell_dataset}
with a row for every unique record in circumcision survey data for a given
area. Also includes empirical estimates for circumcision estimates for each
unique record.}

\item{area_lev}{PSNU area level for specific country.}

\item{k_dt_age}{Age knot spacing in spline definitions, Default: 5}

\item{k_dt_time}{Time knot spacing in spline definitions, set to NULL to
disable temporal splines, Default: NULL}

\item{inc_time_tmc}{Indicator variable which decides whether to include
temporal random effects for TMC as well as MMC, Default: FALSE}
}
\value{
List of design matrices for fixed and random effects for medical
and traditional circumcision.
}
\description{
Create design matrices for fixed effects and temporal, age,
spatial and random effects, for both medical and traditional circumcision.
}
\seealso{
\code{\link[threemc]{create_shell_dataset}}
\code{\link[splines]{splineDesign}}
\code{\link[mgcv]{tensor.prod.model.matrix}}
\code{\link[methods]{as}}
\code{\link[Matrix]{sparse.model.matrix}}
}
\keyword{internal}
