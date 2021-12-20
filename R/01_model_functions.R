#### Modelling Functions ####

#' @title Create Precision Matrix for RWp Process
#' 
#' @description Create the precision matrix for a RWp process.
#' 
#' @param dim Dimension of the precision matrix.
#' @param order Order of the random walk.
#' @param offset.diag Option to offset diagonal by 1E-6.
#' 
#' @return RW precision matrix
#' @export
create_rw_prec_matrix <- function(dim = NULL,
                                  order = 1,
                                  offset.diag = TRUE) {
  # Creating stucture matrix
  Q <- diff(diag(dim), differences = order)
  Q <- t(Q) %*% Q
  # Adding offset to diagonal if required
  if (offset.diag) {
    diag(Q) <- diag(Q) + 1E-6
  }
  # Converting to sparse matrix
  Q <- as(Q, "sparseMatrix")
  # Returning matrix
  return(Q)
}


#' @title Create Precision Matrix for ICAR Process
#' 
#' @description Create the precision matrix for an ICAR process.
#' 
#' @param sf_obj Shapefiles needed for adjacency.
#' @param row.names Unique IDs for the areas.
#' 
#' @return ICAR precision matrix
#' @export
create_icar_prec_matrix <- function(sf_obj = NULL,
                                    row.names = NULL) {
  # Creating neighbourhood structure
  Q_space <-  poly2nb(sf_obj, row.names = sf_obj[, row.names])
  # Converting to adjacency matrix
  Q_space <- nb2mat(Q_space, style = 'B', zero.policy = TRUE)
  # Converting to sparse matrix
  Q_space <- as(Q_space, "sparseMatrix")
  # Creating precision matrix from adjacency
  Q_space <- INLA::inla.scale.model(
    diag(rowSums(Q_space)) - 0.99 * Q_space,
    constr = list(A = matrix(1, 1, nrow(Q_space)), e = 0)
  )
}


#' @title Create Matrix to Estimate Lagged Cumulative Hazard Rate
#'
#' @description Create a matrix to estimate the lagged cumulative hazard rate 
#' needed for survival analysis by age and time. The option to include an 
#' additional stratification variable is also available, creating a 3D hazard 
#' function.
#' 
#' @param dat Dataset used for modelling.
#' 
#' @param time1 Variable name for time of birth.
#' @param subset Subset for dataset.
#' @param time2 Variable name for time circumcised or censored.
#' @param timecaps Window to fix temporal dimension before and after.
#' @param Ntime Number of time points (if NULL, function will calculate).
#' @param age - Variable with age circumcised or censored.
#' @param Nage Number of age groups (if NULL, function will calculate).
#' @param strat Variable to stratify by in using a 3D hazard function.
#' @param Nstrat Number of stratification groups (if NULL, function will 
#' calculate).
#' 
#' @return Matrix for selecting instananeous hazard rate.
#' @export
create_integration_matrix_agetime_lag <- function(dat,
                                                  subset = NULL,
                                                  time1 = 'time1',
                                                  time2 = 'time2',
                                                  timecaps = c(1,Inf),
                                                  Ntime = NULL,
                                                  age = 'age',
                                                  Nage = NULL,
                                                  strat = NULL,
                                                  Nstrat = NULL) {
  # Integration matrix for cumulative hazard
  dat$time1_cap <- pmin(timecaps[2] - timecaps[1] + 1, 
                        pmax(1, as.numeric(dat[[time1]]) - timecaps[1] + 1))
  # Integration matrix for cumulative hazard
  dat$time2_cap <- pmin(timecaps[2] - timecaps[1] + 1, 
                        pmax(1, as.numeric(dat[[time2]]) - timecaps[1] + 1))
  
  # Shifting time points by the time caps
  dat$time1_cap2 <- dat[[time1]] - timecaps[1] + 1
  dat$time2_cap2 <- dat[[time2]] - timecaps[1] + 1
  
  # Number of dimensions in the hazard function
  if (is.null(Ntime)) Ntime <- max(dat[,'time1_cap', drop = TRUE])
  if (is.null(Nage)) Nage <- max(dat[age])
  if (!is.null(strat) & is.null(Nstrat)) Nstrat <- max(dat[strat])
  # Subsetting data if necessary
  if (!is.null(subset)) {
    dat <- subset(dat, eval(parse(text = subset)))
  }
  # Number of rows in the resulting matrix
  nrow <- nrow(dat)
  # Adding dummy variable for the rows of the matrix
  dat$row <- seq_len(nrow(dat))
  
  # Matrix for 3D hazard function if strat not NULL
  if (is.null(strat)) {
    
    # column entries for integration matrix
    cols <- apply(dat, 1, function(x) {
      # If circumcised at birth select relevant entry
      if (as.numeric(x['time1_cap2']) == (as.numeric(x['time2_cap2']))) {
        test <- min(timecaps[2] - timecaps[1] + 1,
                    max(1, as.numeric(x['time1_cap2'])))
      } else {
        # Else just estimate the
        test <- cumsum(c(max(1, as.numeric(x['time1_cap2'])),
                         Ntime + (as.numeric(x['time1_cap2']):(as.numeric(x['time2_cap2']) - 1) > 0 & 
                                    as.numeric(x['time1_cap2']):(as.numeric(x['time2_cap2']) - 1) <= timecaps[2] - timecaps[1])))
      }
      (test <- test[-length(test)])
    })
    cols <- unlist(cols)
    
    # Row entries for integration matrix
    rows <- apply(dat, 1, function(x) {
      rep(as.numeric(x['row']), as.numeric(x[time2]) - as.numeric(x[time1]))
    })
    rows <- unlist(rows)
    
    ncol <- Ntime * Nage
  }
  # Matrix for 3D hazard function if strat not NULL
  if (!is.null(strat)) {
    
    # column entries for integration matrix
    cols <- apply(dat, 1, FUN = function(x) {
      # If circumcised at birth select relevant entry
      if (as.numeric(x['time1_cap2']) == (as.numeric(x['time2_cap2']))) {
        test <- Ntime * Nage * (as.numeric(x[strat]) - 1) +
          min(timecaps[2] - timecaps[1] + 1,
              max(1, as.numeric(x['time1_cap2'])))
      } else {
        # Else just estimate the ?
        test <- cumsum(c(Ntime * Nage * (as.numeric(x[strat]) - 1)
                         + max(1, as.numeric(x['time1_cap2'])),
                         Ntime + (as.numeric(x['time1_cap2']):(as.numeric(x['time2_cap2']) - 1) > 0 & 
                                    as.numeric(x['time1_cap2']):(as.numeric(x['time2_cap2']) - 1) <= timecaps[2] - timecaps[1])))
      }
      test <- test[-length(test)]
      return(test)
    })
    cols <- unlist(cols)
    
    # Row entries for integration matrix
    rows <- apply(dat, 1, function(x) {
      rep(as.numeric(x['row']), as.numeric(x[time2]) - as.numeric(x[time1]))
    })
    rows <- unlist(rows)
    
    ncol <- Ntime * Nage * Nstrat
  }
  # Outputting sparse matrix
  A <- sparseMatrix(i = rows,
                    j = cols,
                    x = 1,
                    dims = c(nrow, ncol))
  # Returning matrix
  return(A)
}
