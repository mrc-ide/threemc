
# function to create integration matrices for selecting instantaneous hazard
# rate

#' @title Create Integration Matrices for Selecting Instantaneous Hazard
#' @description Create (unlagged and lagged) integration matrices for selecting 
#' instantaneous hazard rate, required for fitting TMB model.
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
#' Default: NULL
#' @param  ... Further arguments passed to or from other methods.
#' @return `list` of length 2 of integration matrices for selecting 
#' instantaneous hazard rate.
#' 
#' @seealso
#'  \code{\link[threemc]{create_shell_dataset}}
#'  \code{\link[threemc]{create_integration_matrix_agetime}}
#'  \code{\link[threemc]{create_integration_matrix_agetime_lag}}
#' @rdname create_integration_matrices
#' @export
#' 
#' @importFrom threemc create_integration_matrix_agetime create_integration_matrix_agetime_lag

create_integration_matrices <- function(out,
                                        time1 = "time1",
                                        time2 = "time2",
                                        age = "age",
                                        strat = "space",
                                        ...
                                        ) {
  
  out <- out
  out$time1 <- out$time - out$circ_age
  out$time2 <- out$time
  out$age <- out$circ_age + 1
  
  # Matrix for selecting instantaneous hazard rate
  IntMat1 <- create_integration_matrix_agetime(
    dat = out,
    time1 = time1,
    time2 = time2,
    strat = strat,
    age = age,
    Ntime = length(unique(out$time)),
    ...
  )
  
  # Matrix for selecting instantaneous hazard rate
  IntMat2 <- create_integration_matrix_agetime_lag(
    dat = out,
    time1 = "time1",
    time2 = "time2",
    strat = "space",
    age = "age",
    Ntime = length(unique(out$time)),
    ...
  )
  
  output <- list(
    "IntMat1" = IntMat1,
    "IntMat2" = IntMat2
  )
  return(output)
}
