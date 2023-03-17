
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
#' @param k_dt_age Age knot spacing in spline definitions, Default: 5
#' @param k_dt_time Time knot spacing in spline definitions, set to NULL to 
#' disable temporal splines, Default: NULL
#' @param paed_age_cutoff Age at which to split MMC design matrices between
#' paediatric and non-paediatric populations, the former of which are constant
#' over time. Set to NULL if not desired, Default: NULL
#' @param rw_order Order of the random walk used for temporal precision matrix.
#' Setting to NULL assumes you wish to specify an AR 1 temporal prior.
#' Default: NULL
#' @param inc_time_tmc Indicator variable which decides whether to include
#' temporal random effects for TMC as well as MMC, Default: FALSE
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
#' @importFrom R.utils Arguments
#' @export
threemc_prepare_model_data <- function(out,
                                       areas,
                                       # options
                                       area_lev = NULL,
                                       aggregated = TRUE,
                                       weight = "population",
                                       k_dt_age = 5,
                                       k_dt_time = NULL,
                                       paed_age_cutoff = NULL,
                                       rw_order = NULL,
                                       inc_time_tmc = FALSE,
                                       ...) {
  if (is.null(area_lev)) {
    message(
      "area_lev arg missing, taken as maximum area level in shell dataset"
    )
    area_lev <- max(out$area_level, na.rm = TRUE)
  }
  
  type_info <- TRUE
  if (all(out$obs_mmc == 0) && all(out$obs_tmc == 0)) {
    message("No circumcision type information present in `out`")
    type_info <- FALSE
  }

  # Create design matrices for fixed effects and temporal, age, space and
  # interaction random effects
  design_matrices <- create_design_matrices(
    dat          = out, 
    area_lev     = area_lev, 
    k_dt_age     = k_dt_age, 
    k_dt_time    = k_dt_time, 
    inc_time_tmc = inc_time_tmc
  )

  # Have piecewise mmc design matrices for paediatric and non-paediatric pops
  if (!is.null(paed_age_cutoff)) {
    design_matrices <- split_mmc_design_matrices_paed(
      out, area_lev, design_matrices, paed_age_cutoff
    )
  }

  # Create integration matrices for selecting the instantaneous hazard rate
  integration_matrices <- create_integration_matrices(
    out,
    area_lev = area_lev,
    time1    = "time1",
    time2    = "time2",
    age      = "age",
    strat    = "space",
    ...
  )

  # create survival matrices for MMC, TMC, censored and left censored
  survival_matrices <- create_survival_matrices(
    out,
    areas      = areas,
    area_lev   = area_lev,
    time1      = "time1",
    time2      = "time2",
    age        = "age",
    strat      = "space",
    aggregated = aggregated,
    weight     = weight,
    ...
  )

  # Precision/Adjacency matrix for the spatial random effects
  Q_space <- list(
    "Q_space" = create_icar_prec_matrix(
      sf_obj = areas, area_lev = area_lev, row.names = "space"
    )
  )

  # returned list
  dat_tmb <- c(
    design_matrices, integration_matrices, survival_matrices, Q_space
  )

  # Precision matrix for temporal random effects
  if (!is.null(rw_order)) {
    stopifnot(rw_order %in% c(1, 2))
    message("Random Walk ", rw_order, " temporal prior specified")
    Q_time <- list(
      "Q_time" = create_rw_prec_matrix(
        dim = ncol(design_matrices$X_time_mmc), # same dims as Q_space
        order = rw_order,
        ...
      )
    )
    dat_tmb <- c(dat_tmb, Q_time)
  } else {
    message("rw_order = NULL, AR 1 temporal prior specified")
  }

  # Combine Data for TMB model (also add whether type info is present in out)
  return(c(dat_tmb, "type_info" = type_info))
}

#### create_design_matrices ####

#' @title Create Design Matrices
#' @description Create design matrices for fixed effects and temporal, age,
#' spatial and random effects, for both medical and traditional circumcision.
#' @inheritParams threemc_prepare_model_data
#' @param dat Shell dataset (datputted by \link[threemc]{create_shell_dataset}
#' with a row for every unique record in circumcision survey data for a given
#' area. Also includes empirical estimates for circumcision estimates for each
#' unique record.
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
create_design_matrices <- function(dat,
                                   area_lev = NULL,
                                   k_dt_age = 5,
                                   k_dt_time = 5,
                                   inc_time_tmc = FALSE) {
  if (is.null(area_lev)) {
    message(
      "area_lev arg missing, taken as maximum area level in shell dataset"
    )
    area_lev <- max(dat$area_level, na.rm = TRUE)
  }

  # Only doing the matrices on the specified aggregation
  dat <- shell_data_spec_area(dat, area_lev)

  # Spline definitions
  k_age <- k_dt_age * (floor(min(dat$age) / k_dt_age) - 3):
    (ceiling(max(dat$age) / k_dt_age) + 3)
  if (!is.null(k_dt_time)) {
    k_time <- k_dt_time * (floor(min(dat$time) / k_dt_time) - 3):
      (ceiling(max(dat$time) / k_dt_time) + 3)
  }

  # Design matrix for the fixed effects
  X_fixed <- Matrix::sparse.model.matrix(N ~ 1, data = dat)

  # Design matrix for the temporal random effects
  if (is.null(k_dt_time)) {
    X_time <- Matrix::sparse.model.matrix(N ~ -1 + as.factor(time), data = dat)
  } else {
    X_time <- splines::splineDesign(k_time, dat$time, outer.ok = TRUE)
    X_time <- methods::as(X_time, "sparseMatrix")
  }

  # Design matrix for the age random effects
  X_age <- splines::splineDesign(k_age, dat$age, outer.ok = TRUE)
  X_age <- methods::as(X_age, "sparseMatrix")

  # Design matrix for the spatial random effects
  if (all(dat$space == 1)) {
    form <- stats::formula(N ~ 1)
  } else {
    form <- stats::formula(N ~ -1 + as.factor(space))
  }
  X_space <- Matrix::sparse.model.matrix(form, data = dat)

  # Design matrix for the interaction random effects
  X_agetime <- mgcv::tensor.prod.model.matrix(list(X_time, X_age))
  X_agespace <- mgcv::tensor.prod.model.matrix(list(X_space, X_age))
  if (is.null(k_dt_time)) {
    space_time_fct <- dat %>%
      dplyr::group_by(space, time) %>%
      dplyr::group_indices() %>%
      as.factor()
    X_spacetime <- Matrix::sparse.model.matrix(
      dat$N ~ -1 + space_time_fct
    )
  } else {
    X_spacetime <- mgcv::tensor.prod.model.matrix(list(X_space, X_time))
  }
  
  # return design matrices as list
  output <- list(
    "X_fixed_mmc"     = X_fixed,
    "X_fixed_tmc"     = X_fixed,
    "X_time_mmc"      = X_time,
    "X_time_tmc"      = X_time,
    "X_age_mmc"       = X_age,
    "X_age_tmc"       = X_age,
    "X_space_mmc"     = X_space,
    "X_space_tmc"     = X_space,
    "X_agetime_mmc"   = X_agetime,
    "X_agespace_mmc"  = X_agespace,
    "X_agespace_tmc"  = X_agespace,
    "X_spacetime_mmc" = X_spacetime
  )

  if (inc_time_tmc == FALSE) output <- output[names(output) != "X_time_tmc"]

  return(output)
}

#### split_mmc_design_matrices_paed ####

#' @title Split MMC Design Matrices between adult and paediatric populations
#'
#' @description Take MMC-related design matrices and effectively "split"
#' them into two parts; one paediatric set of design matrices which does not
#' include any time-related components, and one non-paediatric set.
#'
#' @inheritParams threemc_prepare_model_data
#' @param design_matrices Design matrices for fixed effects and temporal, age,
#' spatial and random effects, for both medical and traditional circumcision.
#' @importFrom rlang .data
#' @rdname split_mmc_design_matrices_paed
#' @keywords internal
split_mmc_design_matrices_paed <- function(out,
                                           area_lev,
                                           design_matrices,
                                           paed_age_cutoff = 10) {

  # TODO: What order should these design matrices be in?
  if (!is.null(paed_age_cutoff)) {
    # pull out for area level of interest
    out_spec_area_lev <- dplyr::filter(out, .data$area_level == area_lev)
    # identify rows corresponding to paediatric and adult circumcisions
    paed_age_rows <- which(out_spec_area_lev$circ_age < paed_age_cutoff)
    adult_age_rows <- which(out_spec_area_lev$circ_age >= paed_age_cutoff)

    # pull mmc-related design matrices
    design_matrices_mmc <- design_matrices[
      grepl("mmc", names(design_matrices))
    ]
    # pull non-time matrices; these will be paediatric mmc design matrices
    design_matrices_mmc_paed <- design_matrices_mmc[
      !grepl("time", names(design_matrices_mmc))
    ]
    # append names with "_paed"
    names(design_matrices_mmc_paed) <- paste0(
      names(design_matrices_mmc_paed), "_paed"
    )

    # in adult mmc design matrices, set all rows for paediatric circs to 0
    design_matrices_mmc_adult <- lapply(design_matrices_mmc, function(x) {
      x[paed_age_rows, ] <- 0
      return(x)
    })
    # do the opposite for paediatric mmc design matrices
    design_matrices_mmc_paed <- lapply(design_matrices_mmc_paed, function(x) {
      x[adult_age_rows, ] <- 0
      return(x)
    })

    # append mmc design matrices with adult and paediatric mmc matrices
    design_matrices[grepl("mmc", names(design_matrices))] <-
      design_matrices_mmc_adult
    design_matrices <- c(design_matrices, design_matrices_mmc_paed)
  }
  return(design_matrices)
}


#### shell_data_spec_area ####

#' @title Subset and Prepare Shell Dataset to Only One Admin Boundary
#'
#' @description  Subset the shell dataset to a specific area of interest and
#' resets the space counter. The output dataset can be used to set up model
#' components on the specified administrative boundaries.
#'
#' @param dat Shell dataset used for modelling
#' @inheritParams threemc_prepare_model_data
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname shell_data_spec_area
#' @keywords internal
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
    # Reset counter on space
    dplyr::mutate(space = .data$space - min(.data$space) + 1)

  # Return matrix
  return(dat)
}


#### create_integration_matrices ####

#' @title Create Integration Matrices for Selecting Instantaneous Hazard
#' @description Create (unlagged and lagged) integration matrices for selecting
#' instantaneous hazard rate, required for fitting TMB model.
#' @inheritParams threemc_prepare_model_data
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

  # Prepare age and time variables
  out$time1 <- out$time - out$circ_age
  out$time2 <- out$time
  out$age <- out$circ_age + 1

  # Matrix for selecting instantaneous hazard rate
  # Note: must be kept in CamelCase to agree with syntax for tmb::MakeAdFun
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


#### create_integration_matrix_agetime ####

#' @title Create Matrix to Estimate Cumulative Hazard Rate
#'
#' @description Create a matrix to estimate the cumulative hazard rate needed
#' for survival analysis by age and time. The option to include an additional
#' stratification variable is also available, creating a 3D hazard function.
#'
#' @inheritParams create_integration_matrices
#' @inheritParams create_integration_matrix_agetime
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

  # Shift time points by the time caps
  dat$time1_cap2 <- dat[[time1]] - timecaps[1] + 1
  dat$time2_cap2 <- dat[[time2]] - timecaps[1] + 1

  # If no stratification variable create a dummy variable
  if (is.null(strat)) {
    strat <- "strat"
    dat$strat <- 1
  }

  # Number of dimensions in the hazard function
  if (is.null(Ntime)) Ntime <- max(dat[["time1_cap"]])
  if (is.null(Nage)) Nage <- max(dat[age])
  if (is.null(strat) == FALSE && is.null(Nstrat)) Nstrat <- max(dat[strat])

  # Subset data if necessary
  if (is.null(subset) == FALSE) {
    dat <- subset(dat, eval(parse(text = subset)))
  }
  # Add dummy variable for the rows of the matrix
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

  # Output sparse matrix
  A <- Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = 1,
    dims = c(nrow(dat), ncol)
  )
  # Return matrix
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
#' @inheritParams create_integration_matrices
#' @inheritParams create_integration_matrix_agetime
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

  # Shift time points by the time caps
  dat$time1_cap2 <- dat[[time1]] - timecaps[1] + 1
  dat$time2_cap2 <- dat[[time2]] - timecaps[1] + 1

  # If no stratification variable create a dummy variable
  if (is.null(strat)) {
    strat <- "strat"
    dat$strat <- 1
  }

  # Number of dimensions in the hazard function
  if (is.null(Ntime)) Ntime <- max(dat[, "time1_cap", drop = TRUE])
  if (is.null(Nage)) Nage <- max(dat[age])
  if (!is.null(strat) && is.null(Nstrat)) Nstrat <- max(dat[strat])
  # Subset data if necessary
  if (!is.null(subset)) {
    dat <- subset(dat, eval(parse(text = subset)))
  }
  # Number of rows in the resulting matrix
  nrow <- nrow(dat)

  # Add dummy variable for the rows of the matrix
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

  # Output sparse matrix
  A <- Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = 1,
    dims = c(nrow, ncol)
  )
  # Return matrix
  return(A)
}

#### create_survival_matrices ####

#' @title Create Survival Matrices for Selecting Instantaneous Hazard
#' @description Create empirical agetime hazard matrices for medical and
#' traditional circumcision, as well as for censored (i.e. non-circumcised) and
#' left-censored (i.e. age of circumcision unknown) individuals.
#'
#' @inheritParams threemc_prepare_model_data
#' @inheritParams create_integration_matrices
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

  # calculate empirical agetime hazard matrices for different circ types
  circs <- c(
    "obs_mmc", # medical circumcision rate
    "obs_tmc", # traditional circumcision rate,
    "obs_mc", # all circumcision (to model unknown type)
    "cens", # censored
    "icens" # left censored
  )
  # survival matrix names in TMB model
  list_names <- c("A_mmc", "A_tmc", "A_mc", "B", "C")

  # don't model for specific type if missing from out
  # if all out[[circs[i]]] are 0, create dummy survival matrix of all 0s
  is_missing <- is_dummy <- c()
  dummy_hazard_matrices <- vector(mode = "list", length = length(circs))
  for (i in seq_along(circs)) {
    if (!circs[i] %chin% names(out)) {
      message(paste0("Not creating survival matrix for type == ", circs[i]))
      is_missing <- c(is_missing, i)
    } else if (all(out[[circs[[i]]]] == 0)) {
      message(paste0("Producing dummy survival matrix for type == ", circs[i]))
      dummy_hazard_matrices[[i]] <- Matrix::sparseMatrix(
        i = 1,
        j = nrow(out),
        x = 0,
        dims = c(1, nrow(out))
      )
      is_dummy <- c(is_dummy, i)
    }
  }
  if (!is.null(is_missing) || !is.null(is_dummy)) {
    circs <- circs[-c(is_missing, is_dummy)]
    dummy_names <- list_names[is_dummy]
    if (!is.null(is_missing)) list_names <- list_names[-is_missing]
  }
  # remove any NULL dummy hazard matrices
  dummy_hazard_matrices <- dummy_hazard_matrices[-which(
    vapply(dummy_hazard_matrices, is.null, logical(1))
  )]

  # Matrices for selecting instantaneous hazard rate for:
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

  if (length(dummy_hazard_matrices) != 0) {
    for (i in seq_along(dummy_hazard_matrices)) {
      hazard_matrices <- append(
        hazard_matrices, dummy_hazard_matrices[[i]], is_dummy[i]
      )
    }
  }
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
#' @inheritParams threemc_prepare_model_data 
#' @inheritParams create_design_matrices
#' @inheritParams create_integration_matrix_agetime
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

  # Integration matrix for cumulative hazard
  dat$time2_cap <- pmin(
    timecaps[2] - timecaps[1] + 1,
    pmax(1, as.numeric(dat[[time2]]) - timecaps[1] + 1)
  )

  # If no stratification variable create a dummy variable
  if (is.null(strat)) {
    strat <- "strat"
    dat$strat <- 1
  }

  # Number of dimensions in the hazard function
  if (is.null(Ntime)) Ntime <- max(dat[, "time1_cap", drop = TRUE])
  if (is.null(Nage)) Nage <- max(dat[age])
  if (is.null(Nstrat)) Nstrat <- max(dat[strat])

  # Subsett data if necessary
  if (!is.null(subset)) {
    dat <- subset(dat, eval(parse(text = subset)))
  }

  # If the selection matrices need to be taken from one reference aggregation
  # then we get a list of the hierarchical structure to that level
  if (aggregated == TRUE) {
    # If no weight variable create a dummy variable
    if (is.null(weight)) {
      weight <- "weight"
      dat$weight <- 1
    }

    # Get aggregation structure
    areas_agg <- create_aggregate_structure(
      areas = areas,
      area_lev = area_lev
    )

    # Merge on number of times to
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

    # Only keep strata where we have data
    dat2 <- subset(dat, eval(parse(text = paste(circ, " != 0", sep = "")))) %>%
      dplyr::mutate(row = seq_len(dplyr::n()))

    # Aggregation for each row in the dataframe
    entries <- apply(dat2, 1, function(x) {
      # Get areas in reference administrative
      # boundaries to aggregate over
      tmp_space <- areas_agg$sub_region_list[[as.numeric(x[strat])]]
      # Get columns with non-zero entries for sparse matrix
      cols <- Ntime * Nage * (tmp_space - min_ref_space) +
        Ntime * (as.numeric(x[age]) - 1) +
        as.numeric(x["time2_cap"])
      # Get rows for sparse matrix
      rows <- rep(as.numeric(x["row"]), length(cols))
      # Get weights
      vals <-
        as.numeric(x[circ]) *
          dat[cols, weight, drop = TRUE] / sum(dat[cols, weight, drop = TRUE])
      # Output dataset
      tmp <- data.frame(cols, rows, vals)
      # Return dataframe
      return(tmp)
    })

    # Extract entries for sparse matrix
    cols <- as.numeric(unlist(lapply(entries, "[", "cols")))
    rows <- as.numeric(unlist(lapply(entries, "[", "rows")))
    vals <- as.numeric(unlist(lapply(entries, "[", "vals")))

    # Else the selection matrices will be taken from the aggregation they are on
  } else {
    # Only keep strata where we have data
    dat2 <- subset(dat, eval(parse(text = paste(circ, " != 0", sep = ""))))

    # Column entries for hazard matrix
    cols <- unlist(apply(dat2, 1, function(x) {
      Ntime * Nage * (as.numeric(x[strat]) - 1) + Ntime *
        (as.numeric(x[age]) - 1) + as.numeric(x["time2_cap"])
    }, simplify = FALSE))
    rows <- seq_len(nrow(dat2))
    vals <- dat2[[circ]]
  }

  # sparse haz matrix which selects corresponding incidence rates for lik.
  A <- Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = vals,
    dims = c(nrow(dat2), Ntime * Nage * Nstrat)
  )

  # Return matrix
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
    # Create neighbourhood structure
    Q_space <- spdep::poly2nb(sf_obj, row.names = sf_obj[, row.names])
    # Convert to adjacency matrix
    Q_space <- spdep::nb2mat(Q_space, style = "B", zero.policy = TRUE)

    # for precision matrix
    Q <- diag(rowSums(as.matrix(Q_space))) - 0.99 * Q_space
  } else {
    Q_space <- Matrix::Matrix(data = 0, nrow = 1, ncol = 1)
    Q <- as.matrix(0)
  }

  # Convert to sparse matrix
  Q_space <- methods::as(Q_space, "sparseMatrix")

  # Create precision matrix from adjacency
  Q_space <- scale_gmrf_precision(
    Q   = Q,
    A   = matrix(1, 1, nrow(Q_space)),
    eps = 0
  )

  # Change to same class as outputted by INLA::inla.scale.model
  Q_space <- methods::as(Q_space, "dgTMatrix")
}

#' @title Scale of GMRF precision matrix
#' @description See also `naomi::scale_gmrf_precision()`.
#' @rdname scale_gmrf_precision
#' @keywords internal
scale_gmrf_precision <- function(
    Q, 
    A = matrix(1, ncol = ncol(Q)), 
    eps = sqrt(.Machine$double.eps)
  ) {
  nb <- spdep::mat2listw(abs(Q))$neighbours
  comp <- spdep::n.comp.nb(nb)
  for (k in seq_len(comp$nc)) {
    idx <- which(comp$comp.id == k)
    Qc <- Q[idx, idx, drop = FALSE]
    if (length(idx) == 1) {
      Qc[1, 1] <- 1
    } else {
      qinv <- function(Q, A = NULL) {
        Sigma <- Matrix::solve(Q)
        if (is.null(A)) {
          return(Sigma)
        } else {
          A <- matrix(1, 1, nrow(Sigma))
          W <- Sigma %*% t(A)
          Sigma_const <- Sigma - W %*% Matrix::solve(A %*% W) %*% Matrix::t(W)
          return(Sigma_const)
        }
      }
      
      Ac <- A[, idx, drop = FALSE]
      Qc_eps <- Qc + Matrix::Diagonal(ncol(Qc)) * max(Matrix::diag(Qc)) * 
        eps
      Qc_inv <- qinv(Qc_eps, A = Ac)
      scaling_factor <- exp(mean(log(Matrix::diag(Qc_inv))))
      Qc <- scaling_factor * Qc
    }
    Q[idx, idx] <- Qc
  }
  Q
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
#' @rdname create_rw_prec_matrix
#' @keywords internal
create_rw_prec_matrix <- function(dim,
                                  order = 1,
                                  offset.diag = TRUE) {
  # Create structure matrix
  Q <- diff(diag(dim), differences = order)
  Q <- t(Q) %*% Q
  # Add offset to diagonal if required
  if (offset.diag) {
    diag(Q) <- diag(Q) + 1E-6
  }
  # Convert to sparse matrix
  Q <- methods::as(Q, "sparseMatrix")
  # Return matrix
  return(Q)
}
