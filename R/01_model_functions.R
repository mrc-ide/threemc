#' @title Function to read in Circumcision Data
#' 
#' @description Function to read in circumcision data to fit model. Handles 
#' csvs with `data.table::fread`, and geographical data with `sf::read_sf` (for
#' which it also adds unique identifiers for each `area_level`).
#' 
#' @param path Path to data.
#' @param filters Optional named vector, whose values dictate the values 
#' filtered for in the corresponding column names. Only supports filtering for 
#' one value for each column.
#' 
#' @return relevant data set, filtered as desired.
#' @export
#' @import dplyr
#' @importFrom data.table fread
#' @importFrom sf read_sf
#' @importFrom rlang sym
# imports data.table, sf, dplyr, rlang
# maybe add a warning for missing "circ" columns for surveys?? And add them in
# NAs in this situation (look at KEN for this)
read_circ_data <- function(path, filters = NULL) {
  
  # read in data, depending on file type
  if (grepl(".geojson", path)) {
    .data <- read_sf(path)
  } else .data <- as.data.frame(fread(path))
  
  # if desired, recursively filter data with provided `filters` vector
  if (!is.null(filters)) {
    cols <- names(filters)
    vals <- as.vector(filters[seq_along(filters)])
    for (i in seq_along(filters)) {
      if (!cols[i] %in% names(.data)) next
      .data <- .data %>% 
        # change col i to symbol (if present), evaluate corresponding filter
        dplyr::filter({{ sym(cols[i]) }} == vals[i])
    }
  }
  # for areas, add unique identifier within Admin code and merge to boundaries
  if ("sf" %in% class(.data)) {
    .data <- .data %>%
      group_by(area_level) %>%
      mutate(space = row_number()) %>%
      ungroup()
  }
  return(.data)
}


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
create.rw.prec.matrix <- function(dim = NULL,
                                  order = 1,
                                  offset.diag = TRUE) {
  # Creating stucture matrix
  Q <- diff(diag(dim), differences = order)
  Q <- t(Q) %*% Q
  # Adding offset to diagonal if required
  if (offset.diag == TRUE) {
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
create.icar.prec.matrix <- function(sf_obj = NULL,
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


#' @title Create Matrix Selecting Instananeous Hazard Rate
#' 
#' @description Create a matrix selecting the instantaneous hazard rate needed 
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
#' @param circ Variables with circumcision matrix.
#' 
#' @return Matrix for selecting instananeous hazard rate.
#' @export
create.hazard.matrix.agetime <- function(dat,
                                         subset = NULL,
                                         time1 = 'time1',
                                         time2 = 'time2',
                                         timecaps = c(1, Inf),
                                         Ntime = NULL,
                                         age = 'age',
                                         Nage = NULL,
                                         strat = NULL,
                                         Nstrat = NULL,
                                         circ = 'circ') {
  
  # Integration matrix for cumululative hazard
  dat$time1_cap <- pmin(timecaps[2] - timecaps[1] + 1, 
                        pmax(1, as.numeric(dat[[time1]]) - timecaps[1] + 1))
  
  # Integration matrix for cumululative hazard
  dat$time2_cap <- pmin(timecaps[2] - timecaps[1] + 1, 
                        pmax(1, as.numeric(dat[[time2]]) - timecaps[1] + 1))
  
  # Number of dimensions in the hazard function
  if (is.null(Ntime)) Ntime <- max(dat[,'time1_cap', drop = TRUE])
  if (is.null(Nage)) Nage <- max(dat[age])
  if (is.null(strat) == FALSE & is.null(Nstrat)) Nstrat <- max(dat[strat])
  
  # Subsetting data if necessary
  if (is.null(subset) == FALSE) {
    dat <- subset(dat, eval(parse(text = subset)))
  }
  # Matrix for 2D age time hazard function if strat is NULL
  if (is.null(strat) == TRUE) {
    
    cols <- apply(dat, 1, function(x) {
      Ntime * (as.numeric(x[age]) - 1) + as.numeric(x['time2_cap'])
    })
    cols <- unlist(cols)
    
    # Matrix dimension
    ncol <- Ntime * Nage
  }
  # Matrix for 3D hazard function if strat not NULL
  if (is.null(strat) == FALSE) {
    
    # Integration matrix for cumululative hazard
    cols <- apply(dat, 1, function(x) {
      Ntime * Nage * (as.numeric(x[strat]) - 1) + Ntime * 
        (as.numeric(x[age]) - 1) + as.numeric(x['time2_cap'])
    })
    cols <- unlist(cols)
    
    ncol <- Ntime * Nage * Nstrat
  }
  # Outputting sparse matrix
  A <- sparseMatrix(i = 1:nrow(dat),
                    j = cols,
                    x = dat[[circ]],
                    dims = c(nrow(dat), ncol))
  # Returning matrix
  return(A)
}


#'@title Create Matrix to Estimate Cumulative Hazard Rate
#'
#'@description Create a matrix to estimate the cumulative hazard rate needed
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
create.integration.matrix.agetime <- function(dat,
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
  
  # Integration matrix for cumululative hazard
  dat$time1_cap <- pmin(timecaps[2] - timecaps[1], 
                        pmax(1, as.numeric(dat[[time1]]) - timecaps[1] + 1))
  
  # Integration matrix for cumululative hazard
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
  dat$row <- 1:nrow(dat)
  
  # Matrix for 3D hazard function if strat not NULL
  if (is.null(strat) == TRUE) {
    
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
  if (is.null(strat) == FALSE) {
    
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


#'@title Create Matrix to Estimate Lagged Cumulative Hazard Rate
#'
#'@description Create a matrix to estimate the lagged cumulative hazard rate 
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
#' @param age - Variable with age circumcisied or censored.
#' @param Nage Number of age groups (if NULL, function will calculate).
#' @param strat Variable to stratify by in using a 3D hazard function.
#' @param Nstrat Number of stratification groups (if NULL, function will 
#' calculate).
#' 
#' @return Matrix for selecting instananeous hazard rate.
#' @export
create.integration.matrix.agetime.lag <- function(dat,
                                                  subset = NULL,
                                                  time1 = 'time1',
                                                  time2 = 'time2',
                                                  timecaps = c(1,Inf),
                                                  Ntime = NULL,
                                                  age = 'age',
                                                  Nage = NULL,
                                                  strat = NULL,
                                                  Nstrat = NULL) {
  # Integration matrix for cumululative hazard
  dat$time1_cap <- pmin(timecaps[2] - timecaps[1] + 1, 
                        pmax(1, as.numeric(dat[[time1]]) - timecaps[1] + 1))
  # Integration matrix for cumululative hazard
  dat$time2_cap <- pmin(timecaps[2] - timecaps[1] + 1, 
                        pmax(1, as.numeric(dat[[time2]]) - timecaps[1] + 1))
  
  # Shifting time points by teh time caps
  dat$time1_cap2 <- dat[[time1]] - timecaps[1] + 1
  dat$time2_cap2 <- dat[[time2]] - timecaps[1] + 1
  
  # Number of dimensions in the hazard function
  if (is.null(Ntime)) Ntime <- max(dat[,'time1_cap', drop = TRUE])
  if (is.null(Nage)) Nage <- max(dat[age])
  if (is.null(strat) == FALSE & is.null(Nstrat)) Nstrat <- max(dat[strat])
  # Subsetting data if necessary
  if (is.null(subset) == FALSE) {
    dat <- subset(dat, eval(parse(text=subset)))
  }
  # Nnumber of rows in the resulting matrix
  nrow <- nrow(dat)
  # Adding dummy variabel for the rows of the matrix
  dat$row <- 1:nrow(dat)
  # Matrix for 3D hazard function if strat not NULL
  if (is.null(strat) == TRUE) {
    
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
  if (is.null(strat) == FALSE) {
    
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


#' @title Prepare Survey Data
#' 
#' @description Prepare survey data required to run the circumcision model.
#' 
#' @param area_hierarchy - Hierarchy and metadata of administrative boundaries.
#' @param area_boundaries - SF file with administrative boundaries (shapefiles).
#' @param survey_circumcision - Information on male circumcision status from 
#' surveys.
#' @param survey_individuals - Information on the individuals surveyed.
#' @param survey_clusters - Information on the survey clusters.
#' @param area_lev - Desired admin boundary level to perform the analysis on.
#' @param start_year - Year to begin the analysis on.
#' @param cens_year - Year to censor the circumcision data by (Sometimes some 
#' weirdness at the final survey year, e.g. v small number of MCs).
#' @param cens_age - Age to censor the circumcision data at (Default 59, 
#' i.e. no circumcisions after 59 years of age).
#' 
#' @return Survey data with required variables to run circumcision model.
#' @export
prepare.survey.data <- function(areas,
                                survey_circumcision,
                                survey_individuals,
                                survey_clusters,
                                area_lev = 0,
                                start_year = 2006,
                                cens_year = NULL,
                                cens_age = 59) {
  
  ## Merging circumcision and individuals survey datasets ---------------------
  
  # Bringing datasets together
  survey_circumcision <- survey_circumcision %>%
    # Merging on individual information to  the circumcision dataset
    dplyr::left_join(
      survey_individuals %>%
        dplyr::select(survey_id, cluster_id, individual_id, sex, age, indweight),
      by = c("survey_id", "individual_id")) %>%
    # Merging on cluster information to the circumcision dataset
    dplyr::left_join(
      survey_clusters %>%
        dplyr::select(c(survey_id, cluster_id, area_id = geoloc_area_id)),
      by = c("survey_id", "cluster_id")) %>%
    # Remove those with missing circumcison status
    dplyr::filter(!is.na(circ_status) & !is.na(age) & !is.na(indweight)) %>%
    # Variables needed for analysis
    dplyr::mutate(
      # Survey year
      year = as.numeric(substr(survey_id, 4, 7)),
      # Year of Birth (estimated as no DOB filly yet)
      yob = year - age,
      # If circumcision age > age of the individual set, reset circumcision age
      circ_age = ifelse(circ_age > age, NA, circ_age)
    )
  
  ##############################
  ## Censoring if necessary ###
  ##############################
  # Censoring at cens_year if assumed no circumcisions after a certain year
  if (is.null(cens_age) == FALSE) {
    survey_circumcision <- survey_circumcision %>%
      # Censoring indivduals from analysis at cens_age
      mutate(
        # No circumcision after cens_age
        circ_status = ifelse(circ_status == 1 & !is.na(circ_age) & 
                               circ_age > cens_age, 0.0, circ_status),
        # Resetting age at circumcision
        circ_age = ifelse(circ_age > cens_age, NA, circ_age),
        # Resetting age for everyone else
        age = ifelse(age > cens_age, cens_age, age),
        # Year of circ/censoring (estimated using the age as no date of circ)
        yoc = ifelse(!is.na(circ_age), yob + circ_age, yob + age)
      )
  }
  
  # Censoring at cens_year if assumed no circumcisions after a certain year
  if (is.null(cens_year) == FALSE) {
    survey_circumcision <- survey_circumcision %>%
      # Censoring at cens_year
      dplyr::filter(yob < cens_year) %>%
      # Final variables for modelling
      dplyr::mutate(
        # Censoring circumcision status for those circumcised in cens_year,
        # Assuming the interval censored people were circumcised before cens_year
        circ_status = ifelse(yoc >= cens_year & circ_status == 1 & 
                               !is.na(circ_age), 0.0, circ_status),
        # circ censoring year (or censor year in cens_year - 1) at cens_year - 1
        yoc = ifelse(yoc == cens_year, cens_year - 1, yoc)
      )
  }
  
  #########################################
  ### Setting desired level aggregation ###
  #########################################
  # Getting the area level id to province
  for (i in seq_len(max(areas$area_level))) {
    survey_circumcision <- survey_circumcision %>%
      # Merging on boundary information
      dplyr::left_join(
        areas %>%
          sf::st_drop_geometry() %>%
          dplyr::select(area_id, area_level, parent_area_id),
        by = "area_id"
      ) %>%
      # Altering area
      mutate(area_id = if_else(area_level == area_lev,
                               as.character(area_id),
                               as.character(parent_area_id))) %>%
      dplyr::select(-c(parent_area_id, area_level))
  }
  
  ###################################################
  ### Final preparation of circumcision variables ###
  ###################################################
  # Preparing circumcision variables for the model
  survey_circumcision <- survey_circumcision %>%
    # Merging on the region index
    dplyr::left_join(
      areas %>%
        sf::st_drop_geometry() %>%
        dplyr::select(area_id, area_name, space),
      by = "area_id"
    ) %>%
    dplyr::mutate(
      # Time interval for the individual
      time1 = yob - start_year + 1,
      time2 = yoc - start_year + 1,
      # Event type
      event = ifelse(circ_status == 1 & !is.na(circ_age), 1,
                     ifelse((circ_status == 1 & is.na(circ_age)), 2, 0)),
      # Circumcision age
      circ_age = yoc - yob,
      age = circ_age + 1)
  
  # Returning prepped circumcision datasets
  return(survey_circumcision)
}


#' @title Normalise Survey Weights and apply Kish Coefficients
#' 
#' @description Normalise survey weights and apply Kish coefficients.
#' 
#' @param survey_circumcision Information on male circumcision status from
#' surveys containing survey weights.
#' @param strata.norm Stratification variables for normalising survey weights.
#' @param strara.kish Stratification vairables for estimating and applying the 
#' Kish coefficients.
#'  
#' @return Survey data with normalised survey weights and required variables to 
#' run circumcision model.
#' @export
normalise.weights.kish <- function(survey_circumcision,
                                   strata.norm = c('survey_id', 'area_id'),
                                   strata.kish = c('survey_id')) {
  # Preparing survey weights for the model
  survey_circumcision <- survey_circumcision %>%
    # Standardising survey weights
    group_by_at(.vars = strata.norm) %>%
    mutate(indweight_st = indweight / mean(indweight, na.rm=TRUE)) %>%
    ungroup() %>%
    # Applying Kish coefficient to the survey weights
    left_join(
      plyr::ddply(survey_circumcision,
                  .variables = strata.kish,
                  summarize,
                  N = length(survey_id),
                  Neff = (sum(indweight) ^2)/sum(indweight * indweight),
                  ratio = N/Neff),
      by = 'survey_id') %>%
    mutate(indweight_st = indweight_st / ratio)
}
