#' @title Create Survival Matrices for Selecting Instantaneous Hazard
#' @description Create empirical agetime hazard matrices for medical and
#' traditional circumcision, as well as for censored (i.e. non-circumcised) and
#' left-censored (i.e. age of circumcision unknown) individuals.
#'
#' @param out Shell dataset (outputted by \link[threemc]{create_shell_dataset}
#' with a row for every unique record in circumcision survey data for a given
#' area. Also includes empirical estimates for circumcision estimates for each
#' unique record.
#' @param areas `sf` shapefiles for specific country/region.
#' @param area_lev  PSNU area level for specific country. Defaults to the
#' maximum area level found in `areas` if not supplied.
#' @param time1 Variable name for time of birth, Default: "time1"
#' @param time2 Variable name for time circumcised or censored,
#' Default: "time2"
#' @param age - Variable with age circumcised or censored. Default: "age"
#' @param strat Variable to stratify by in using a 3D hazard function,
#' Default: "space"
#' @param aggregated ??
#' @param  ... Further arguments passed to or from other methods.
#' @return `list` of length 4 of survival matrices for selecting
#' instantaneous hazard rate.
#'
#' @seealso
#'  \code{\link[threemc]{create_shell_dataset}}
#'  \code{\link[threemc]{create_hazard_matrix_agetime}}
#' @rdname create_survival_matrices
#' @export
create_survival_matrices <- function(out,
                                     areas = areas,
                                     area_lev = area_lev,
                                     time1 = "time1",
                                     time2 = "time2",
                                     age = "age",
                                     strat = "space",
                                     aggregated = TRUE,
                                     ...) {
  out$time1 <- out$time - out$circ_age
  out$time2 <- out$time
  
  ## calculate empirical agetime hazard matrices for different circ types
  circs <- c(
    "obs_mmc", # medical circumcision rate
    "obs_tmc", # traditional circumcision rate,
    "obs_mc", # all circumcision (to model unknown type)
    "cens", # censored
    "icens" # left censored
  )
  # suffices for entries in the final list 
  list_names <- c("_mmc", "_tmc", "_mc", "_rc", "_lc")
  # remove MC if modelling for missing type is undesirable
  if (!"obs_mc" %in% names(out)) {
    circs <- circs[-3]
    list_names <- list_names[-3]
  }
  ## Matrices for selecting instantaneous hazard rate for:
  hazard_matrices <- lapply(circs, function(x) {
    threemc::create_hazard_matrix_agetime(
      dat = out,
      areas = areas,
      area_lev = area_lev,
      time1 = time1,
      time2 = time2,
      strat = strat,
      age   = age,
      circ  = x,
      Ntime = length(unique(out$time)),
      aggregated = TRUE,
      ...
    )
  })
  # Turning list of five lists into one list 
  hazard_matrices <- rlang::flatten(hazard_matrices)
  # Altering entry names
  names(hazard_matrices) <- paste(names(hazard_matrices), rep(list_names, each = 2), sep = '')
  # Retunring hazard matrices
  return(hazard_matrices)
}
