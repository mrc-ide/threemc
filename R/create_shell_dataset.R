#' @title Create Shell Dataset for Estimating Empirical Circumcision Rate
#'
#' @description  Create a shell dataset with a row for every unique area ID,
#' area name, year and circumcision age in survey data. Also, computes the
#' empirical number of person years until circumcision and number of people
#' circumcised for several "types" of circumcision; known medical
#' circumcisions, known traditional circumcisions, censored survey entries
#' (i.e. where surveyed individuals had not been circumcised) and left-censored
#'  survey entries (i.e. where circumcision occurred at an unknown age).
#'
#' @param survey_circumcision Information on male circumcision status from
#' surveys.
#' @param areas `sf` shapefiles for specific country/region.
#' @param area_lev  PSNU area level for specific country. Defaults to the
#' maximum area level found in `areas` if not supplied.
#' @param time1 Variable name for time of birth, Default: "time1"
#' @param time2 Variable name for time circumcised or censored,
#' Default: "time2"
#' @param strat Variable to stratify by in using a 3D hazard function.
#' @param age - Variable with age circumcised or censored. Default: "age"
#' @param circ Variables with circumcision matrix, Default: "indweight_st"
#' @param ...  Further arguments passed to or from other methods.
#'
#' @seealso
#'  \code{\link[threemc]{datapack_psnu_area_level}}
#'  \code{\link[tidyr]{crossing}}
#'  \code{\link[threemc]{create_integration_matrix_agetime}}
#'  \code{\link[threemc]{create_hazard_matrix_agetime}}
#' @return `data.frame` with a row for every unique record in
#' `survey_circumcision` for a given area. Also includes empirical estimates
#' for circumcision estimates for each unique record.
#' @export
#'
#' @import dplyr
#' @import sf
#' @import rlang
create_shell_dataset <- function(survey_circumcision,
                                 areas,
                                 area_lev,
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

  if (missing(area_lev)) {
    message("area_lev arg missing, taken as maximum area level in areas")
    area_lev <- max(areas$area_level, na.rm = TRUE)
  }

  if (!"Matrix" %in% .packages()) {
      message(paste("Strongly recommend loading 'Matrix' package, as it is",
                    "is required for summing sparse matrices"))
  }


  ## remove spatial elements from areas, take only specified/highest area level
  if (inherits(areas, "sf")) {
    areas_model <- st_drop_geometry(areas)
  } else {
    areas_model <- areas
  }

  areas_model <- areas_model %>%
      filter(.data$area_level == area_lev) %>%
      select(any_of(c("area_id", "area_name", "space")))

  ## create skeleton dataset with row for every unique area_id, area_name,
  ## space, year and circ_age
  out <- tidyr::crossing(areas_model,
    "year" = seq(2006, 2021, by = 1),
    "circ_age" = 0 : max(survey_circumcision$circ_age, na.rm = TRUE)
  ) %>%
    ## Getting time and age variable
    mutate(
      time = .data$year - 2006 + 1,
      age = .data$circ_age + 1
    ) %>%
    ## Sorting dataset
    arrange(.data$space, .data$age, .data$time)

  ## Obtain N person years
  out_int_mat <- threemc::create_integration_matrix_agetime(
    dat = survey_circumcision,
    time1 = time1,
    time2 = time2,
    strat = strat,
    age   = age,
    Ntime = length(unique(out$time)),
    ...
  )
  # browser()
  out$N <- as.vector(survey_circumcision$indweight_st %*% out_int_mat)

  ## calculate empirical agetime hazard matrices for different circumcision
  ## types, and take column sums (i.e. N empirical circs for each "type):
  subsets <- c(
    "event == 1 & type == 'MMC'", # N MMC
    "event == 1 & type == 'TMC'", # N TMC
    "event == 0", # N censored (i.e. not circumcised)
    "event == 2" # N left-censored (circ at unknown age)
  )
  agetime_hazard_matrices <- lapply(subsets, function(x) {
    threemc::create_hazard_matrix_agetime(
      dat = survey_circumcision,
      subset = x,
      time1 = time1,
      time2 = time2,
      strat = strat,
      age   = age,
      circ  = circ,
      Ntime = length(unique(out$time)),
      ...
    ) |>
    colSums()
  })

  ## add to out:
  empirical_circ_cols <- c("obs_mmc", "obs_tmc", "cens", "icens")
  out[, empirical_circ_cols] <- agetime_hazard_matrices

  return(out)
}
