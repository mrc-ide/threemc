
#### Main Function ####

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

#### create_design_matrices #### 

#' @title Create Design Matrices
#' @description Create design matrices for fixed effects and temporal, age,
#' spatial and random effects, for both medical and traditional circumcision.
#'
#' @param dat Shell dataset (datputted by \link[threemc]{create_shell_dataset}
#' with a row for every unique record in circumcision survey data for a given
#' area. Also includes empirical estimates for circumcision estimates for each
#' unique record.
#' @param area_lev Desired admin boundary level to perform the analysis on.
#' @param k_dt Age knot spacing in spline definitions, Default: 5
#' @return List of design matrices for fixed and random effects for medical
#' and traditional circumcision.
#'
#' @seealso
#'  \code{\link[threemc]{create_shell_dataset}}
#'  \code{\link[splines]{splineDesign}}
#'  \code{\link[mgcv]{tensor.prod.model.matrix}}
#'  \code{\link[methods]{as}}
#'  \code{\link[Matrix]{sparse.model.matrix}}
#' @rdname create_design_matrices
#' @importFrom dplyr %>%
#' @keywords internal
create_design_matrices <- function(dat, area_lev = NULL, k_dt = 5) {
  if (is.null(area_lev)) {
    message(
      "area_lev arg missing, taken as maximum area level in shell dataset"
    )
    area_lev <- max(dat$area_level, na.rm = TRUE)
  }

  # Only doing the matrices on the specified aggregation
  dat <- shell_data_spec_area(dat, area_lev)

  ## Spline definitions
  k_age <-
    k_dt * (floor(min(dat$age) / k_dt) - 3):(ceiling(max(dat$age) / k_dt) + 3)

  ## Design matrix for the fixed effects
  X_fixed <- Matrix::sparse.model.matrix(N ~ 1, data = dat)

  ## Design matrix for the temporal random effects
  X_time <- Matrix::sparse.model.matrix(N ~ -1 + as.factor(time), data = dat)

  ## Design matrix for the age random effects
  X_age <- splines::splineDesign(k_age, dat$age, outer.ok = TRUE)
  X_age <- methods::as(X_age, "sparseMatrix")

  ## Design matrix for the spatial random effects
  if (all(dat$space == 1)) {
    form <- stats::formula(N ~ 1)
  } else {
    form <- stats::formula(N ~ -1 + as.factor(space))
  }
  X_space <- Matrix::sparse.model.matrix(form, data = dat)
  
  ## Design matrix for the interaction random effects
  X_agetime <- mgcv::tensor.prod.model.matrix(list(X_time, X_age))
  X_agespace <- mgcv::tensor.prod.model.matrix(list(X_space, X_age))
  X_spacetime <- Matrix::sparse.model.matrix(
    N ~ -1 + factor((dat %>%
      group_by(space, time) %>%
      group_indices())),
    data = dat
  )
  ## return design matrices as list
  output <- list(
    "X_fixed_mmc"     = X_fixed,
    "X_fixed_tmc"     = X_fixed,
    "X_time_mmc"      = X_time,
    "X_age_mmc"       = X_age,
    "X_age_tmc"       = X_age,
    "X_space_mmc"     = X_space,
    "X_space_tmc"     = X_space,
    "X_agetime_mmc"   = X_agetime,
    "X_agespace_mmc"  = X_agespace,
    "X_agespace_tmc"  = X_agespace,
    "X_spacetime_mmc" = X_spacetime
  )
  return(output)
}

#### shell_data_spec_area ####

#' @title Subset and Prepare Shell Dataset to Only One Admin Boundary
#'
#' @description  Subset the shell dataset to a specific area of interest and
#' resets the space counter. The output dataset can be used to set up model
#' components on the specified administrative boundaries.
#'
#' @param dat Shell dataset used for modelling
#' @param area_lev  PSNU area level for specific country. Defaults to the
#' Defaults to the maximum area level found in `dat` if not supplied.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname shell_data_spec_area
shell_data_spec_area <- function(dat, area_lev = NULL) {
  if (is.null(area_lev)) {
    message(
      "area_lev arg missing, taken as maximum area level in shell dataset"
    )
    area_lev <- max(dat$area_level, na.rm = TRUE)
  }

  # Only doing the matrices on the specified aggregation
  dat <- dat %>%
    dplyr::filter(.data$area_level == area_lev) %>%
    # Resetting counter on space
    dplyr::mutate(space = .data$space - min(.data$space) + 1)

  ## Returning matrix
  return(dat)
}


#### create_integration_matrices #### 


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
#' @param area_lev  PSNU area level for specific country. Defaults to the
#' Defaults to the maximum area level found in `dat` if not supplied.
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
#' @keywords internal
create_integration_matrices <- function(out,
                                        area_lev = NULL,
                                        time1 = "time1",
                                        time2 = "time2",
                                        age = "age",
                                        strat = "space",
                                        ...) {
  if (is.null(area_lev)) {
    message(
      "area_lev arg missing, taken as maximum area level in shell dataset"
    )
    area_lev <- max(out$area_level, na.rm = TRUE)
  }

  # Only doing the matrices on the specified aggregation
  out <- shell_data_spec_area(out, area_lev)

  # Preparing age and time variables
  out$time1 <- out$time - out$circ_age
  out$time2 <- out$time
  out$age <- out$circ_age + 1

  ## Matrix for selecting instantaneous hazard rate
  ## Note: must be kept in CamelCase to agree with syntax for tmb::MakeAdFun
  IntMat1 <- create_integration_matrix_agetime(
    dat = out,
    time1 = time1,
    time2 = time2,
    strat = strat,
    age = age,
    Ntime = length(unique(out$time)),
    ...
  )

  ## Matrix for selecting instantaneous hazard rate
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


#### create_integration_matrix_agetime ####

#' @title Create Matrix to Estimate Cumulative Hazard Rate
#'
#' @description Create a matrix to estimate the cumulative hazard rate needed
#' for survival analysis by age and time. The option to include an additional
#' stratification variable is also available, creating a 3D hazard function.
#'
#' @param dat Dataset used for modelling.
#' @param subset Subset for dataset, Default: NULL
#' @param time1 Variable name for time of birth, Default: "time1"
#' @param time2 Variable name for time circumcised or censored, Default: "time2"
#' @param timecaps Window to fix temporal dimension before and after,
#' Default: c(1, Inf)
#' @param Ntime Number of time points (if NULL, function will calculate),
#' Default: NULL
#' @param age - Variable with age circumcised or censored, Default: "age"
#' @param Nage Number of age groups (if NULL, function will calculate),
#' Default: NULL
#' @param strat Variable to stratify by in using a 3D hazard function,
#' Default: NULL
#' @param Nstrat Number of stratification groups (if NULL, function will
#' calculate), Default: NULL
#
#' @return Matrix for selecting instantaneous hazard rate.
#' @seealso
#'   \code{\link[Matrix]{sparseMatrix}}
#' @rdname create_integration_matrix_agetime
#' @keywords internal
create_integration_matrix_agetime <- function(dat,
                                              subset = NULL,
                                              time1 = "time1",
                                              time2 = "time2",
                                              timecaps = c(1, Inf),
                                              Ntime = NULL,
                                              age = "age",
                                              Nage = NULL,
                                              strat = NULL,
                                              Nstrat = NULL) {


  # !! JE: Matt -- check these lines; I think this can be done with
  # pmin()/pmax() and does not need an unlist() because it will always return
  # a vector.

  # Integration matrix for cumulative hazard
  dat$time1_cap <- pmin(
    timecaps[2] - timecaps[1],
    pmax(1, as.numeric(dat[[time1]]) - timecaps[1] + 1)
  )

  # Integration matrix for cumulative hazard
  dat$time2_cap <- pmin(
    timecaps[2] - timecaps[1] + 1,
    pmax(1, as.numeric(dat[[time2]]) - timecaps[1] + 1)
  )

  # Shifting time points by the time caps
  dat$time1_cap2 <- dat[[time1]] - timecaps[1] + 1
  dat$time2_cap2 <- dat[[time2]] - timecaps[1] + 1

  ## If no stratification variable create a dummy variable
  if (is.null(strat)) {
    strat <- "strat"
    dat$strat <- 1
  }

  # Number of dimensions in the hazard function
  if (is.null(Ntime)) Ntime <- max(dat[["time1_cap"]])
  if (is.null(Nage)) Nage <- max(dat[age])
  if (is.null(strat) == FALSE & is.null(Nstrat)) Nstrat <- max(dat[strat])

  # Subsetting data if necessary
  if (is.null(subset) == FALSE) {
    dat <- subset(dat, eval(parse(text = subset)))
  }
  # Adding dummy variable for the rows of the matrix
  dat$row <- seq_len(nrow(dat))

  # column entries for integration matrix
  cols <- unlist(apply(dat, 1, function(x) {
    # If circumcised at birth select relevant entry
    if (as.numeric(x["time1_cap2"]) == (as.numeric(x["time2_cap2"]))) {
      Ntime * Nage * (as.numeric(x[strat]) - 1) +
        min(
          timecaps[2] - timecaps[1] + 1,
          max(1, as.numeric(x["time1_cap2"]))
        )
    } else {
      # Else just estimate the
      cumsum(
        c(
          Ntime * Nage * (as.numeric(x[strat]) - 1) +
            max(1, as.numeric(x["time1_cap2"])),
          Ntime + (as.numeric(x["time1_cap2"]):
          (as.numeric(x["time2_cap2"]) - 1) > 0 &
            as.numeric(x["time1_cap2"]):
            (as.numeric(x["time2_cap2"]) - 1) <=
              timecaps[2] - timecaps[1])
        )
      )
    }
  }, simplify = FALSE))

  # Matrix dimension
  ncol <- Ntime * Nage * Nstrat

  # Row entries for integration matrix
  rows <- unlist(apply(dat, 1, function(x) {
    rep(as.numeric(x["row"]), as.numeric(x[time2]) - as.numeric(x[time1]) + 1)
  }, simplify = FALSE))

  # Outputting sparse matrix
  A <- Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = 1,
    dims = c(nrow(dat), ncol)
  )
  # Returning matrix
  return(A)
}

#### create_integration_matrix_agetime_lag #### 

#' @title Create Matrix to Estimate Lagged Cumulative Hazard Rate
#'
#' @description Create a matrix to estimate the lagged cumulative hazard rate
#' needed for survival analysis by age and time. The option to include an
#' additional stratification variable is also available, creating a 3D hazard
#' function.
#'
#' @param dat Dataset used for modelling.
#' @param subset Subset for dataset, Default: NULL
#' @param time1 Variable name for time of birth, Default: "time1"
#' @param time2 Variable name for time circumcised or censored, Default: "time2"
#' @param timecaps Window to fix temporal dimension before and after,
#' Default: c(1, Inf)
#' @param Ntime Number of time points (if NULL, function will calculate),
#' Default: NULL
#' @param age - Variable with age circumcised or censored, Default: "age"
#' @param Nage Number of age groups (if NULL, function will calculate),
#' Default: NULL
#' @param strat Variable to stratify by in using a 3D hazard function,
#' Default: NULL
#' @param Nstrat Number of stratification groups (if NULL, function will
#' calculate), Default: NULL
#' @return Matrix for selecting instantaneous hazard rate.
#'
#' @seealso
#'  \code{\link[Matrix]{sparseMatrix}}
#' @rdname create_integration_matrix_agetime_lag
#' @keywords internal
create_integration_matrix_agetime_lag <- function(dat,
                                                  subset = NULL,
                                                  time1 = "time1",
                                                  time2 = "time2",
                                                  timecaps = c(1, Inf),
                                                  Ntime = NULL,
                                                  age = "age",
                                                  Nage = NULL,
                                                  strat = NULL,
                                                  Nstrat = NULL) {
  # Integration matrix for cumulative hazard
  dat$time1_cap <- pmin(
    timecaps[2] - timecaps[1] + 1,
    pmax(1, as.numeric(dat[[time1]]) - timecaps[1] + 1)
  )
  # Integration matrix for cumulative hazard
  dat$time2_cap <- pmin(
    timecaps[2] - timecaps[1] + 1,
    pmax(1, as.numeric(dat[[time2]]) - timecaps[1] + 1)
  )

  # Shifting time points by the time caps
  dat$time1_cap2 <- dat[[time1]] - timecaps[1] + 1
  dat$time2_cap2 <- dat[[time2]] - timecaps[1] + 1

  ## If no stratification variable create a dummy variable
  if (is.null(strat)) {
    strat <- "strat"
    dat$strat <- 1
  }

  # Number of dimensions in the hazard function
  if (is.null(Ntime)) Ntime <- max(dat[, "time1_cap", drop = TRUE])
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

  # column entries for integration matrix
  cols <- unlist(apply(dat, 1, FUN = function(x) {
    # If circumcised at birth select relevant entry
    if (as.numeric(x["time1_cap2"]) == (as.numeric(x["time2_cap2"]))) {
      test <- Ntime * Nage * (as.numeric(x[strat]) - 1) +
        min(
          timecaps[2] - timecaps[1] + 1,
          max(1, as.numeric(x["time1_cap2"]))
        )
    } else {
      # Else just estimate the
      test <- cumsum(
        c(
          Ntime * Nage * (as.numeric(x[strat]) - 1) +
            max(1, as.numeric(x["time1_cap2"])),
          Ntime + (as.numeric(x["time1_cap2"]):
          (as.numeric(x["time2_cap2"]) - 1) > 0 &
            as.numeric(x["time1_cap2"]):
            (as.numeric(x["time2_cap2"]) - 1) <=
              timecaps[2] - timecaps[1])
        )
      )
    }
    test <- test[-length(test)]
    return(test)
  }, simplify = FALSE))

  # Number of columns
  ncol <- Ntime * Nage * Nstrat

  # Row entries for integration matrix
  rows <- unlist(apply(dat, 1, function(x) {
    rep(as.numeric(x["row"]), as.numeric(x[time2]) - as.numeric(x[time1]))
  }, simplify = FALSE))

  # Outputting sparse matrix
  A <- Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = 1,
    dims = c(nrow, ncol)
  )
  # Returning matrix
  return(A)
}

#### create_survival_matrices ####

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
#' @param aggregated `agggregated = FALSE` treats every area_id as its own
#' object, allowing for the use of surveys for lower area hierarchies. 
#' `aggregated = TRUE` means we only look at area level of interest.
#' @param weight variable to weigh circumcisions by when aggregating for 
#' lower area hierarchies (only applicable for `aggregated = TRUE`) 
#' @param  ... Further arguments passed to or from other methods.
#' @return `list` of length 4 of survival matrices for selecting
#' instantaneous hazard rate.
#'
#' @seealso
#'  \code{\link[threemc]{create_shell_dataset}}
#'  \code{\link[threemc]{create_hazard_matrix_agetime}}
#' @rdname create_survival_matrices
#' @keywords internal
create_survival_matrices <- function(out,
                                     areas = areas,
                                     area_lev = area_lev,
                                     time1 = "time1",
                                     time2 = "time2",
                                     age = "age",
                                     strat = "space",
                                     aggregated = TRUE,
                                     weight = "population",
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
  list_names <- c("A_mmc", "A_tmc", "A_mc", "B", "C")
  # remove MC if modelling for missing type is undesirable
  if (!"obs_mc" %in% names(out)) {
    circs <- circs[-3]
    list_names <- list_names[-3]
  }
  ## Matrices for selecting instantaneous hazard rate for:
  hazard_matrices <- lapply(circs, function(x) {
    create_hazard_matrix_agetime(
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
      weight = weight,
      ...
    )
  })
  names(hazard_matrices) <- list_names

  return(hazard_matrices)
}

#### create_hazard_matrix_agetime #### 

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
#' @param aggregated `agggregated = FALSE` treats every area_id as its own
#' object, allowing for the use of surveys for lower area hierarchies. 
#' `aggregated = TRUE` means we only look at area level of interest.
#' @param weight variable to weigh circumcisions by when aggregating for 
#' lower area hierarchies (only applicable for `aggregated = TRUE`) 
#' @return Matrix for selecting instantaneous hazard rate.
#'
#' @seealso
#'  \code{\link[threemc]{create_shell_dataset}}
#' @rdname create_hazard_matrix_agetime
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @keywords internal
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

  # Integration matrix for cumulative hazard
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
  if (is.null(strat)) {
    strat <- "strat"
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
  if (aggregated == TRUE) {
    ## If no weighting variable create a dummy variable
    if (is.null(weight)) {
      weight <- "weight"
      dat$weight <- 1
    }

    # Getting aggregation structure
    areas_agg <- create_aggregate_structure(
      areas = areas,
      area_lev = area_lev
    )

    # Merging on number of times to
    # replicate to the main dataset
    dat <- dat %>%
      dplyr::left_join(
        areas_agg$n_sub_region_df,
        by = "area_id"
      )

    # Minimum space ID within the reference level
    min_ref_space <- min(dat %>%
      dplyr::filter(.data$area_level == area_lev) %>%
      dplyr::pull(.data$space))

    # Minimum space ID within the reference level
    Nstrat <- dat %>%
      dplyr::filter(.data$area_level == area_lev) %>%
      dplyr::distinct(.data$space) %>%
      dplyr::pull() %>%
      length()

    # Only keeping strata where we have data
    dat2 <- subset(dat, eval(parse(text = paste(circ, " != 0", sep = "")))) %>%
      dplyr::mutate(row = 1:dplyr::n())

    # Aggregation for each row in the dataframe
    entries <- apply(dat2, 1, function(x) {
      # Getting areas in reference administrative
      # boundaries to aggregate over
      tmp_space <- areas_agg$sub_region_list[[as.numeric(x[strat])]]
      # Getting columns with non-zero entries for sparse matrix
      cols <- Ntime * Nage * (tmp_space - min_ref_space) +
        Ntime * (as.numeric(x[age]) - 1) +
        as.numeric(x["time2_cap"])
      # Getting rows for sparse matrix
      rows <- rep(as.numeric(x["row"]), length(cols))
      # Getting weights
      vals <-
        as.numeric(x[circ]) *
          dat[cols, weight, drop = TRUE] / sum(dat[cols, weight, drop = TRUE])
      # Output dataset
      tmp <- data.frame(cols, rows, vals)
      # Return dataframe
      return(tmp)
    })

    # Extracting entries for sparse matrix
    cols <- as.numeric(unlist(lapply(entries, "[", "cols")))
    rows <- as.numeric(unlist(lapply(entries, "[", "rows")))
    vals <- as.numeric(unlist(lapply(entries, "[", "vals")))
    
  # Else the selection matrices will be taken from the aggregation they are on
  } else {
    # Only keeping strata where we have data
    dat2 <- subset(dat, eval(parse(text = paste(circ, " != 0", sep = ""))))

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

#### create_icar_prec_matrix ####

#' @title Create Precision Matrix for ICAR Process
#'
#' @description Create the precision matrix for an ICAR process.
#'
#' @param sf_obj Shapefiles needed for adjacency, Default: NULL
#' @param area_lev  PSNU area level for specific country.
#' @param row.names Unique IDs for the areas, Default: NULL
#' @return ICAR precision matrix.
#'
#' @seealso
#'  \code{\link[spdep]{poly2nb}}
#'  \code{\link[spdep]{nb2mat}}
#'  \code{\link[naomi]{scale_gmrf_precision}}
#'  \code{\link[methods]{as}}
#' @rdname create_icar_prec_matrix
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @keywords internal
create_icar_prec_matrix <- function(sf_obj = NULL,
                                    area_lev = NULL,
                                    row.names = NULL) {
  if (is.null(area_lev)) {
    message(
      "area_lev missing, taken as maximum area level in sf_obj"
    )
    area_lev <- max(sf_obj$area_level, na.rm = TRUE)
  }
 
  sf_obj <- sf_obj %>%
    dplyr::filter(.data$area_level == area_lev)
  
  # if area_lev == 0, adjacency matrix will be a 1x1 matrix with single entry 0
  if (area_lev > 0) {
    # Creating neighbourhood structure
    Q_space <- spdep::poly2nb(sf_obj, row.names = sf_obj[, row.names])
    # Converting to adjacency matrix
    Q_space <- spdep::nb2mat(Q_space, style = "B", zero.policy = TRUE)
    
    # for precision matrix
    Q <- diag(rowSums(as.matrix(Q_space))) - 0.99 * Q_space
  } else {
    Q_space <- Matrix::Matrix(data = 0, nrow = 1, ncol = 1)
    Q <- as.matrix(0)
  }
  
  # Converting to sparse matrix
  Q_space <- methods::as(Q_space, "sparseMatrix")
  
  # Creating precision matrix from adjacency
  Q_space <- naomi::scale_gmrf_precision(
    Q   = Q, 
    A   = matrix(1, 1, nrow(Q_space)),
    eps = 0
  )
  
  # Change to same class as outputted by INLA::inla.scale.model
  Q_space <- methods::as(Q_space, "dgTMatrix")
}

#### create_rw_prec_matrix ####

#' @title Create Precision Matrix for RWp Process
#'
#' @description Create the precision matrix for a RWp process.
#'
#' @param dim Dimension of the precision matrix.
#' @param order Order of the random walk, Default: 1
#' @param offset.diag Option to offset diagonal by 1E-6, Default: TRUE
#'
#' @seealso
#'  \code{\link[methods]{as}}
#
#' @return RW precision matrix
#' @keywords internal
create_rw_prec_matrix <- function(dim,
                                  order = 1,
                                  offset.diag = TRUE) {
  ## Creating structure matrix
  Q <- diff(diag(dim), differences = order)
  Q <- t(Q) %*% Q
  ## Adding offset to diagonal if required
  if (offset.diag) {
    diag(Q) <- diag(Q) + 1E-6
  }
  ## Converting to sparse matrix
  Q <- methods::as(Q, "sparseMatrix")
  ## Returning matrix
  return(Q)
}
