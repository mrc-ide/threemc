#' @title Create Survival Matrices for Selecting Instantaneous Hazard 
#' @description Create empirical agetime hazard matrices for medical and 
#' traditional circumcision, as well as for censored (i.e. non-circumcised) and
#' left-censored (i.e. age of circumcision unknown) individuals.
#'
#' @param out Shell dataset (outputted by \link[threemc]{create_shell_dataset}
#' with a row for every unique record in circumcision survey data for a given 
#' area. Also includes empirical estimates for circumcision estimates for each 
#' unique record.
#' @param time1 Variable name for time of birth, Default: "time1"
#' @param time2 Variable name for time circumcised or censored,
#' Default: "time2"
#' @param age - Variable with age circumcised or censored. Default: "age"
#' @param strat Variable to stratify by in using a 3D hazard function,
#' Default: "space"
#' @param  ... Further arguments passed to or from other methods.
#' @return `list` of length 4 of survival matrices for selecting 
#' instantaneous hazard rate.
#' 
#' @seealso
#'  \code{\link[threemc]{create_shell_dataset}}
#'  \code{\link[threemc]{create_hazard_matrix_agetime}}
#' @rdname create_survival_matrices
#' @export
create_survival_matrices <- function(
  out, 
  time1 = "time1",
  time2 = "time2",
  age = "age",
  strat = "space",
  ...) {
  
  out$time1 <- out$time - out$circ_age
  out$time2 <- out$time
  
  ## calculate empirical agetime hazard matrices for different circ types
  circs <- c(
    "obs_mmc", # medical circumcision rate
    "obs_tmc", # traditional circumcision rate
    "cens",    # censored
    "icens"    # left censored
  )
  hazard_matrices <- lapply(circs, function(x) {
    threemc::create_hazard_matrix_agetime(
      dat = out,
      time1 = time1,
      time2 = time2,
      strat = strat,
      age   = age,
      circ  = x,
      Ntime = length(unique(out$time)),
      ...
    )
  })
  
  ## Matrices for selecting instantaneous hazard rate for:
  output <- list(
    "A_mmc" = hazard_matrices[[1]], # MMC 
    "A_tmc" = hazard_matrices[[2]], # TMC
    "B"     = hazard_matrices[[3]], # censored
    "C"     =  hazard_matrices[[4]] # left-censored
  )
  return(output)
}
