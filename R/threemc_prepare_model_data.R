#' @title Produce Data Matrices for Modelling
#' @description Create data for modelling. Output detailed below.
#' @param out Shell dataset (outputted by \link[threemc]{create_shell_dataset}
#' with a row for every unique record in circumcision survey data for a given
#' area. Also includes empirical estimates for circumcision estimates for each
#' unique record.
#' @param areas `sf` shapefiles for specific country/region.
#' @param area_lev PSNU area level for specific country. 
#' @param aggregated `agggregated = FALSE` treats every area_id as its own
#' object, allowing for the use of surveys for lower area hierarchies. 
#' `aggregated = TRUE` means we only look at area level of interest.
#' @param weight variable to weigh circumcisions by when aggregating for 
#' lower area hierarchies (only applicable for `aggregated = TRUE`) 
#' @param k_dt Age knot spacing in spline definitions, Default: 5
#' @param ... Additional arguments to be passed to functions which create
#' matrices. 
#' @return \code{list} of data required for model fitting, including:
#' \itemize{
#'   \item{design_matrices}{Includes`X_fixed_mmc`, `X_fixed_tmc`, `X_time_mmc`, 
#'   `X_age_mmc`, `X_age_tmc`, `X_space_mmc`, `X_space_tmc`, `X_agetime_mmc`,
#'   `X_agespace_mmc`, `X_agespace_tmc`, `X_spacetime_mmc`. Design 
#'    Create design matrices for fixed effects and temporal, age, space and
#'    interaction random effects}
#'    \item{integration matrices}{Includes `IntMat1`, `IntMat2`. Integration 
#'    matrices for selecting the instantaneous hazard rate.}
#'    \item{survival matrices}{Includes `A_mmc`, `A_tmc`, `A_mc`, `B`, `C`. 
#'    Survival matrices for MMC, TMC, censored and left censored}
#'    \item{Q_space}{Precision/Adjacency matrix for the spatial random effects.
#'    }
#' }
#' @rdname threemc_prepare_model_data
#' @export
threemc_prepare_model_data <- function(
  # data
  out, areas, 
  # options
  area_lev, aggregated = TRUE, weight = "population", k_dt = 5, ...
) {
  
  # Create design matrices for fixed effects and temporal, age, space and
  # interaction random effects
  design_matrices <- create_design_matrices(dat = out,
                                            area_lev = area_lev,
                                            k_dt = k_dt)
  
  # Create integration matrices for selecting the instantaneous hazard rate
  integration_matrices <- create_integration_matrices(out,
                                                      area_lev = area_lev,
                                                      time1 = "time1",
                                                      time2 = "time2",
                                                      age = "age",
                                                      strat = "space",
                                                      ...)
  
  # create survival matrices for MMC, TMC, censored and left censored
  survival_matrices <- create_survival_matrices(out,
                                                areas = areas,
                                                area_lev = area_lev,
                                                time1 = "time1",
                                                time2 = "time2",
                                                age = "age",
                                                strat = "space",
                                                aggregated = aggregated,
                                                weight = weight, 
                                                ...)
  
  # Precision/Adjacency matrix for the spatial random effects
  Q_space <- list("Q_space" =
                    create_icar_prec_matrix(sf_obj    = areas,
                                            area_lev  = area_lev,
                                            row.names = "space"))
  
  # Combine Data for tmb model
  return(c(design_matrices, integration_matrices, survival_matrices, Q_space))
}
