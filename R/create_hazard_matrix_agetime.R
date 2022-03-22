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
#' @param areas `sf` shapefiles for specific country/region.
#' @param area_lev  PSNU area level for specific country. 
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
                                         areas, 
                                         area_lev, 
                                         subset = NULL,
                                         time1 = "time1",
                                         time2 = "time2",
                                         timecaps = c(1, Inf),
                                         Ntime = NULL,
                                         age = "age",
                                         Nage = NULL,
                                         strat = NULL,
                                         Nstrat = NULL,
                                         circ = "circ",
                                         aggregated = FALSE,
                                         weight = NULL) {
  
  ## Integration matrix for cumulative hazard
  dat$time1_cap <- pmin(
    timecaps[2] - timecaps[1] + 1,
    pmax(1, as.numeric(dat[[time1]]) - timecaps[1] + 1)
  )
  
  ## Integration matrix for cumulative hazard
  dat$time2_cap <- pmin(
    timecaps[2] - timecaps[1] + 1,
    pmax(1, as.numeric(dat[[time2]]) - timecaps[1] + 1)
  )
  
  ## If no stratification variable create a dummy variable
  if (is.null(strat)){
    strat <- 'strat'
    dat$strat <- 1
  }
  
  ## Number of dimensions in the hazard function
  if (is.null(Ntime)) Ntime <- max(dat[, "time1_cap", drop = TRUE])
  if (is.null(Nage)) Nage <- max(dat[age])
  if (is.null(Nstrat)) Nstrat <- max(dat[strat])
  
  ## Subsetting data if necessary
  if (!is.null(subset)) {
    dat <- subset(dat, eval(parse(text = subset)))
  }
  
  # If the selection matrices need to be taken from one reference aggregation 
  # then we get a list of the hierarchical structure to that level
  if (aggregated == TRUE){
    ## If no weighting variable create a dummy variable
    if (is.null(weight)){
      weight <- 'weight'
      dat$weight <- 1
    }
    
    # Getting aggregation structure
    areas_agg <- create_aggregate_structure(areas    = areas,
                                            area_lev = 5)
    
    # Merging on number of times to 
    # replicate to the main dataset
    dat <- dat %>%
      left_join(areas_agg$areas_agg2,
                by = 'area_id')
    
    # Minimum space ID within the reference level
    min_ref_space <- min(out %>%
                           dplyr::filter(area_level == area_lev) %>%
                           dplyr::pull(space))
    
    # Minimum space ID within the reference level
    Nstrat <- out %>%
      dplyr::filter(area_level == area_lev) %>%
      dplyr::pull(space) %>%
      unique() %>%
      length()
    
    # Only keeping stratums where we have data
    dat2 <- subset(dat, eval(parse(text = paste(circ, ' != 0', sep = '')))) %>%
      dplyr::mutate(row = 1:dplyr::n())
    
    # Aggregation for each row in the dataframe
    entries <- apply(dat2, 1, function(x) {
      # Getting areas in reference administrative 
      # boundaries to aggregate over
      tmp_space <- areas_agg$areas_agg1[[as.numeric(x[strat])]]
      # Getting columns with non-zero entries for sparse matrix
      cols <- Ntime * Nage * (tmp_space - min_ref_space) +
        Ntime * (as.numeric(x[age]) - 1) +
        as.numeric(x["time2_cap"])
      # Getting rows for sparse matrix
      rows <- rep(as.numeric(x["row"]), length(cols))
      # Getting weights 
      vals <- as.numeric(x[circ]) * dat[cols, weight, drop = TRUE] / sum(dat[cols, weight, drop = TRUE])
      # Output dataset 
      tmp <- data.frame(cols, rows, vals)
      # Return dataframe
      return(tmp)
    })
    
    # Extracting entries for sparseMatrix
    cols <- as.numeric(unlist(lapply(entries, "[", "cols")))
    rows <- as.numeric(unlist(lapply(entries, "[", "rows")))
    vals <- as.numeric(unlist(lapply(entries, "[", "vals")))
  }
  # Else the selection matrices will be taken from the aggregation they are on
  else{
    # Only keeping stratums where we have data 
    dat2 <- subset(dat, eval(parse(text = paste(circ, ' != 0', sep = ''))))
    
    ## Column entries for hazard matrix
    cols <- unlist(apply(dat2, 1, function(x) {
      Ntime * Nage * (as.numeric(x[strat]) - 1) + Ntime *
        (as.numeric(x[age]) - 1) + as.numeric(x["time2_cap"])
    }, simplify = FALSE))
    rows <- seq_len(nrow(dat2))
    vals <- dat2[[circ]]
  }
  
  ## Outputting sparse hazard matrix which selects the
  ## corresponding incidence rates for the likelihood.
  A <- Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = vals,
    dims = c(nrow(dat2), Ntime * Nage * Nstrat)
  )
  
  ## Returning matrix
  return(A)
}
