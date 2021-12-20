#' @title Create Matrix to Estimate Cumulative Hazard Rate
#'
#' @description Create a matrix to estimate the cumulative hazard rate needed
#' for survival analysis by age and time. The option to include an additional 
#' stratification variable is also available, creating a 3D hazard function.
#' 
#' @param dat Dataset used for modelling.
#' @param subset Subset for dataset.
#' @param time1 Variable name for time of birth.
#' @param time2 Variable name for time circumcised or censored.
#' @param timecaps Window to fix temporal dimension before and after.
#' @param Ntime Number of time points (if NULL, function will calculate).
#' @param age - Variable with age circumcisied or censored.
#' @param Nage Number of age groups (if NULL, function will calculate).
#' @param strat Variable to stratify by in using a 3D hazard function.
#' @param Nstrat Number of stratification groups (if NULL, function will 
#' calculate).
#' 
#' @return Matrix for selecting instananeous hazard rate.
#' @export
#' 
#' @importFrom Matrix sparseMatrix
create_integration_matrix_agetime <- function(dat,
                                              subset = NULL,
                                              time1 = 'time1',
                                              time2 = 'time2',
                                              timecaps = c(1,Inf),
                                              Ntime = NULL,
                                              age = 'age',
                                              Nage = NULL,
                                              strat = NULL,
                                              Nstrat = NULL) {
  
  
  # !! JE: Matt -- check these lines; I think this can be done with 
  # pmin()/pmax() and does not need an unlist() because it will always return 
  # a vector.
  
  # Integration matrix for cumulative hazard
  dat$time1_cap <- pmin(timecaps[2] - timecaps[1], 
                        pmax(1, as.numeric(dat[[time1]]) - timecaps[1] + 1))
  
  # Integration matrix for cumulative hazard
  dat$time2_cap <- pmin(timecaps[2] - timecaps[1] + 1, 
                        pmax(1, as.numeric(dat[[time2]]) - timecaps[1] + 1))
  
  
  # Shifting time points by the time caps
  dat$time1_cap2 <- dat[[time1]] - timecaps[1] + 1
  dat$time2_cap2 <- dat[[time2]] - timecaps[1] + 1
  
  # Number of dimensions in the hazard function
  if (is.null(Ntime)) Ntime <- max(dat[["time1_cap"]])
  if (is.null(Nage)) Nage <- max(dat[age])
  if (is.null(strat) == FALSE & is.null(Nstrat)) Nstrat <- max(dat[strat])
  
  # Subsetting data if necessary
  if (is.null(subset) == FALSE) {
    dat <- subset(dat, eval(parse(text=subset)))
  }
  # Adding dummy variabld for the rows of the matrix
  dat$row <- seq_len(nrow(dat))
  
  # Matrix for 3D hazard function if strat not NULL
  if (is.null(strat)) {
    
    ## !! JE: Matt -- could this condition be combined with the next one
    ##        by setting Nstrat = 1 if is.null(strat) = TRUE?
    
    # column entries for integration matrix
    cols <- apply(dat, 1, FUN = function(x) {
      # If circumcised at birth select relevant entry
      if (as.numeric(x['time1_cap2']) == (as.numeric(x['time2_cap2']))) {
        min(timecaps[2] - timecaps[1] + 1,
            max(1, as.numeric(x['time1_cap2'])))
        # Else just estimate the ??
      } else {
        cumsum(
          c(max(1, as.numeric(x['time1_cap2'])),
            Ntime + (as.numeric(x['time1_cap2']):(as.numeric(x['time2_cap2']) - 1) > 0 &
                       as.numeric(x['time1_cap2']):(as.numeric(x['time2_cap2']) - 1) <= timecaps[2] - timecaps[1]))
        )
      }
    })
    cols <- unlist(cols)
    
    # Row entries for integration matrix
    rows <- apply(dat, 1, function(x) {
      rep(as.numeric(x['row']), as.numeric(x[time2]) - as.numeric(x[time1]) + 1)
    })
    rows <- unlist(rows)
    
    # Matrix dimension
    ncol <- Ntime * Nage
  }
  # Matrix for 3D hazard function if strat not NULL
  if (!is.null(strat)) {
    
    # column entries for integration matrix
    cols <- apply(dat, 1, function(x) {
      # If circumcised at birth select relevant entry
      if (as.numeric(x['time1_cap2']) == (as.numeric(x['time2_cap2']))) {
        Ntime * Nage * (as.numeric(x[strat]) - 1) +
          min(timecaps[2] - timecaps[1] + 1,
              max(1, as.numeric(x['time1_cap2'])))
      } else {
        # Else just estimate the ?
        cumsum(
          c(Ntime * Nage * (as.numeric(x[strat]) - 1) + max(1, as.numeric(x['time1_cap2'])),
            Ntime + (as.numeric(x['time1_cap2']):(as.numeric(x['time2_cap2']) - 1) > 0 &
                       as.numeric(x['time1_cap2']):(as.numeric(x['time2_cap2']) - 1) <= timecaps[2] - timecaps[1]))
        )
      }
    })
    cols <- unlist(cols)
    
    # Row entries for integration matrix
    rows <- apply(dat, 1, function(x) {
      rep(as.numeric(x['row']), as.numeric(x[time2]) - as.numeric(x[time1]) + 1)
    })
    rows <- unlist(rows)
    
    # Matrix dimension
    ncol <- Ntime * Nage * Nstrat
  }
  
  # Outputting sparse matrix
  A <- sparseMatrix(i = rows,
                    j = cols,
                    x = 1,
                    dims = c(nrow(dat), ncol))
  # Returning matrix
  return(A)
}
