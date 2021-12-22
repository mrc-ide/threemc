#' @title Create Matrix Selecting Instantaneous Hazard Rate
#' 
#' @description Create a matrix selecting the instantaneous hazard rate needed 
#' for survival analysis by age and time. The option to include an additional 
#' stratification variable is also available, creating a 3D hazard function.
#' 
#' @param dat Shell dataset (outputted by \link[threemc]{create_shell_dataset}
#' with a row for every unique record in circumcision survey data for a given 
#' area. Also includes empirical estimates for circumcision estimates for each 
#' unique record.
#' @param subset Subset for dataset, Default: NULL
#' @param time1 Variable name for time of birth, Default: "time1"
#' @param time2 Variable name for time circumcised or censored,
#' Default: "time2"
#' @param timecaps Window to fix temporal dimension before and after,
#' Default: c(1, Inf)
#' @param Ntime Number of time points (if NULL, function will calculate),
#' Default: NULL
#' @param age - Variable with age circumcised or censored. Default: "age"
#' @param Nage Number of age groups (if NULL, function will calculate),
#' Default: NULL
#' @param strat Variable to stratify by in using a 3D hazard function,
#' Default: NULL
#' @param Nstrat Number of stratification groups (if NULL, function will 
#' calculate), Default: NULL
#' @param circ Variables with circumcision matrix, Default: "circ"
#' @return Matrix for selecting instantaneous hazard rate.
#' 
#' @seealso 
#'  \code{\link[threemc]{create_shell_dataset}}
#' @rdname create_hazard_matrix_agetime
#' @export
create_hazard_matrix_agetime <- function(dat,
                                         subset = NULL,
                                         time1 = "time1",
                                         time2 = "time2",
                                         timecaps = c(1, Inf),
                                         Ntime = NULL,
                                         age = "age",
                                         Nage = NULL,
                                         strat = NULL,
                                         Nstrat = NULL,
                                         circ = "circ") {
  
  # Integration matrix for cumulative hazard
  dat$time1_cap <- pmin(timecaps[2] - timecaps[1] + 1, 
                        pmax(1, as.numeric(dat[[time1]]) - timecaps[1] + 1))
  
  # Integration matrix for cumulative hazard
  dat$time2_cap <- pmin(timecaps[2] - timecaps[1] + 1, 
                        pmax(1, as.numeric(dat[[time2]]) - timecaps[1] + 1))
  
  # Number of dimensions in the hazard function
  if (is.null(Ntime)) Ntime <- max(dat[, "time1_cap", drop = TRUE])
  if (is.null(Nage)) Nage <- max(dat[age])
  if (!is.null(strat) & is.null(Nstrat)) Nstrat <- max(dat[strat])
  
  # Subsetting data if necessary
  if (!is.null(subset)) {
    dat <- subset(dat, eval(parse(text = subset)))
  }
  # Matrix for 2D age time hazard function if strat is NULL
  if (is.null(strat)) {
    
    cols <- apply(dat, 1, function(x) {
      Ntime * (as.numeric(x[age]) - 1) + as.numeric(x["time2_cap"])
    })
    cols <- unlist(cols)
    
    # Matrix dimension
    ncol <- Ntime * Nage
  }
  # Matrix for 3D hazard function if strat not NULL
  if (!is.null(strat)) {
    
    # Integration matrix for cumululative hazard
    cols <- apply(dat, 1, function(x) {
      Ntime * Nage * (as.numeric(x[strat]) - 1) + Ntime * 
        (as.numeric(x[age]) - 1) + as.numeric(x["time2_cap"])
    })
    cols <- unlist(cols)
    
    ncol <- Ntime * Nage * Nstrat
  }
  # Outputting sparse matrix
  A <- sparseMatrix(i = seq_len(nrow(dat)),
                    j = cols,
                    x = dat[[circ]],
                    dims = c(nrow(dat), ncol))
  # Returning matrix
  return(A)
}
