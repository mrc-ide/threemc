#' @title Create Shell Dataset for Estimating Empirical Circumcision Rate
#' @description  Create a shell dataset with a row for every unique area ID,
#' area name, year and circumcision age in survey data. Also, computes the
#' empirical number of person years until circumcision and number of people
#' circumcised for several "types" of circumcision; known medical
#' circumcisions, known traditional circumcisions, censored survey entries
#' (i.e. where surveyed individuals had not been circumcised) and left-censored
#'  survey entries (i.e. where circumcision occurred at an unknown age).
#' @param survey_circumcision Information on male circumcision status from
#' surveys.
#' @param populations Single age male population counts by space and time.
#' @param areas `sf` shapefiles for specific country/region.
#' @param area_lev  PSNU area level for specific country. Defaults to the
#' maximum area level found in `areas` if not supplied.
#' @param start_year First year in shell dataset, Default: 2006
#' @param end_year Last year in shell dataset, which is also the year to 
#' forecast/model until, Default: 2021
#' @param time1 Variable name for time of birth, Default: "time1"
#' @param time2 Variable name for time circumcised or censored,
#' Default: "time2"
#' @param strat Variable to stratify by in using a 3D hazard function,
#' Default: "space"
#' @param age - Variable with age circumcised or censored. Default: "age"
#' @param circ Variables with circumcision matrix, Default: "indweight_st"
#' @param ...  Further arguments passed to or from other methods.
#' @seealso
#'  \code{\link[threemc]{datapack_psnu_area_level}}
#'  \code{\link[tidyr]{crossing}}
#'  \code{\link[threemc]{create_integration_matrix_agetime}}
#'  \code{\link[threemc]{create_hazard_matrix_agetime}}
#' @return `data.frame` with a row for every unique record in
#' `survey_circumcision` for a given area. Also includes empirical estimates
#' for circumcision estimates for each unique record.
#' @rdname create_shell_dataset
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr %>%
create_shell_dataset <- function(survey_circumcision,
                                 population_data,
                                 areas,
                                 area_lev = NULL,
                                 start_year = 2006,
                                 end_year = 2021,
                                 time1 = "time1",
                                 time2 = "time2",
                                 strat = "space",
                                 age = "age",
                                 circ = "indweight_st",
                                 ...) {

  ## !! JE: the area_id field here used the area_id that appeared in the
  ##        circumcision dataset. If there were some area_id with no
  ##        observations, then they were dropped, which created a
  ##        misalignment of the indexing.  Instead, use the areas dataset
  ##        to construct this output frame to ensure all districts are
  ##        represented.
  ##
  ##        I suspect that we also want the circ_age field to be constructed
  ##        based on the theoretical maximum circumcision age that we want
  ##        outputs for, rather than the maximum observed age; but not 100%
  ##        sure.

  if (is.null(area_lev)) {
    message("area_lev arg missing, taken as maximum area level in areas")
    area_lev <- max(areas$area_level, na.rm = TRUE)
  }

  ## remove spatial elements from areas, take only specified/highest area level
  if (inherits(areas, "sf")) {
    areas_model <- sf::st_drop_geometry(areas)
  } else {
    areas_model <- areas
  }

  areas_model <- areas_model %>%
    dplyr::filter(.data$area_level <= area_lev) %>%
    dplyr::select(.data$area_id, .data$area_name, .data$area_level, .data$space)

  ## create skeleton dataset with row for every unique area_id, area_name,
  ## space, year and circ_age
  out <- tidyr::crossing(areas_model,
    "year" = seq(start_year, end_year, by = 1),
    "circ_age" = 0:max(survey_circumcision$circ_age, na.rm = TRUE)
  ) %>%
    ## Getting time and age variable
    dplyr::mutate(
      time = .data$year - start_year + 1,
      age = .data$circ_age + 1
    ) %>%
    ## Sorting dataset
    dplyr::arrange(.data$space, .data$age, .data$time) %>%
    ## Adding population data on to merge
    dplyr::left_join(
      populations %>%
        dplyr::select(
          .data$area_id, .data$year,
          circ_age = .data$age, .data$population
        ),
      by = c("area_id", "year", "circ_age")
    )

  ## Add `space` to survey_circumcision observations
  survey_circumcision <- survey_circumcision %>%
    dplyr::left_join(
      dplyr::select(areas_model, .data$area_id, .data$space),
      by = "area_id"
    )
  stopifnot(!is.na(survey_circumcision$space))
  
  ## Obtain N person years
  out_int_mat <- create_integration_matrix_agetime(
    dat = survey_circumcision,
    time1 = time1,
    time2 = time2,
    strat = strat,
    Nstrat = nrow(areas_model),
    age   = age,
    Ntime = length(unique(out$time)),
    ...
  )
  out$N <- as.vector(survey_circumcision$indweight_st %*% out_int_mat)

  ## calculate empirical agetime hazard matrices for different circumcision
  ## types, and take column sums (i.e. N empirical circs for each "type):
  empirical_circ_cols <- c("obs_mmc", "obs_tmc", "obs_mc", "cens", "icens")
  subsets <- c(
    "event == 1 & type == 'MMC'", # N MMC
    "event == 1 & type == 'TMC'", # N TMC
    "event == 1 & type == 'Missing'",
    "event == 0", # N censored (i.e. not circumcised)
    "event == 2" # N left-censored (circ at unknown age)
  )

  # if there is no missing type, there is no need for MC:
  if (!any(survey_circumcision$type == "Missing")) {
    empirical_circ_cols <- empirical_circ_cols[-3]
    subsets <- subsets[-3]
  }

  agetime_hazard_matrices <- lapply(subsets, function(x) {
    create_hazard_matrix_agetime(
      dat = survey_circumcision,
      areas = areas,
      area_lev = area_lev,
      subset = x,
      time1 = time1,
      time2 = time2,
      strat = strat,
      age   = age,
      circ  = circ,
      Ntime = length(unique(out$time)),
      Nstrat = nrow(areas_model),
      ...
    )
  })
  agetime_hazard_matrices <- lapply(agetime_hazard_matrices, Matrix::colSums)

  ## add to out:
  out[, empirical_circ_cols] <- agetime_hazard_matrices
  return(out)
}


#' @title Normalise Survey Weights and apply Kish Coefficients
#' @description Normalise survey weights and apply Kish coefficients.
#' @param survey_circumcision Information on male circumcision status from
#' surveys containing survey weights.
#' @param strata.norm Stratification variables for normalising survey weights,
#' Default: c("survey_id", "area_id")
#' @param strata.kish Stratification variables for estimating and applying the
#' Kish coefficients, Default: "survey_id"
#' @return Survey data with normalised survey weights and required variables to
#' run circumcision model.
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname normalise_weights_kish
#' @keywords internal
normalise_weights_kish <- function(survey_circumcision,
                                   strata.norm = c("survey_id", "area_id"),
                                   strata.kish = c("survey_id")) {

  # Preparing survey weights for the model
  survey_circumcision %>%
    # Standardising survey weights
    dplyr::group_by(dplyr::across(dplyr::all_of(strata.norm))) %>%
    dplyr::mutate(
      indweight_st = .data$indweight / mean(.data$indweight, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    # Applying Kish coefficient to the survey weights
    dplyr::left_join(
      (survey_circumcision %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(strata.kish))) %>%
        dplyr::summarise(
          N = length(.data$survey_id),
          Neff = (sum(.data$indweight)^2) /
            sum(.data$indweight * .data$indweight),
          ratio = .data$N / .data$Neff,
          .groups = "drop"
        )),
      by = "survey_id"
    ) %>%
    dplyr::mutate(indweight_st = .data$indweight_st / .data$ratio)
}
