#' @title Create Shell Dataset for Estimating Empirical Circumcision Rate
#' @description  Create a shell dataset with a row for every unique area ID,
#' area name, year and circumcision age in survey data. Also, computes the
#' empirical number of person years until circumcision and number of people
#' circumcised for several "types" of circumcision; known medical
#' circumcisions, known traditional circumcisions, censored survey entries
#' (i.e. where surveyed individuals had not been circumcised) and left-censored
#'  survey entries (i.e. where circumcision occurred at an unknown age).
#' @inheritParams prepare_survey_data
#' @inheritParams create_integration_matrices
#' @inheritParams threemc_aggregate
#' @param start_year First year in shell dataset.
#' @param end_year Last year in shell dataset, which is also the year to
#' forecast/model until, Default: 2021
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
                                 populations,
                                 areas,
                                 area_lev = NULL,
                                 start_year,
                                 end_year = 2021,
                                 time1 = "time1",
                                 time2 = "time2",
                                 strat = "space",
                                 age = "age",
                                 circ = "indweight_st",
                                 ...) {

  # !! JE: the area_id field here used the area_id that appeared in the
  #        circumcision dataset. If there were some area_id with no
  #        observations, then they were dropped, which created a
  #        misalignment of the indexing.  Instead, use the areas dataset
  #        to construct this output frame to ensure all districts are
  #        represented.
  #
  #        I suspect that we also want the circ_age field to be constructed
  #        based on the theoretical maximum circumcision age that we want
  #        outputs for, rather than the maximum observed age; but not 100%
  #        sure.

  if (is.null(area_lev)) {
    message("area_lev arg missing, taken as maximum area level in areas")
    area_lev <- max(areas$area_level, na.rm = TRUE)
  }

  # check for area_ids with duplicated area_names
  pop_duplicates <- populations %>% 
    dplyr::distinct(.data$area_id, .data$area_name) %>% 
    dplyr::group_by(.data$area_id) %>% 
    dplyr::summarise(n = dplyr::n()) %>% 
    dplyr::arrange(dplyr::desc(.data$n)) %>% 
    filter(.data$n > 1)
  
  if (nrow(pop_duplicates) > 0) {
    
    message("populations has multiple names for the following area_ids:")
    message(paste0(utils::capture.output(
      data.frame(pop_duplicates)), collapse = "\n"
    ))
    message("Removing area_name and treating these as single areas")
    # remove area_name from pops and take unique vals to avoid duplication
    populations <- populations %>% 
      dplyr::select(-.data$area_name) %>% 
      dplyr::distinct()
  }
  
  # check that there is a population for every year
  min_pop_year <- min(populations$year)
  if (start_year < min_pop_year) {
    populations <- fill_downup_populations(
      populations, 
      start_year, 
      min_pop_year
    )
  }
 
  # remove spatial elements from areas, take only specified/highest area level
  if (inherits(areas, "sf")) {
    areas_model <- sf::st_drop_geometry(areas)
  } else {
    areas_model <- areas
  }

  areas_model <- areas_model %>%
    mutate(iso3 = substr(.data$area_id, 0, 3)) %>% 
    dplyr::filter(
      .data$area_level <= area_lev,
      # be sure not to include other countries 
      iso3 %chin% substr(survey_circumcision$area_id, 0, 3)
    ) %>%
    dplyr::select(
      .data$area_id, 
      .data$area_name, 
      contains("parent_area"),
      .data$area_level, 
      .data$space
    )
  
  # "manually" add in parent_area_name, if not in areas
  if (!"parent_area_name" %in% names(areas_model)) {
    parent_areas <- areas_model %>% 
      dplyr::select(
        parent_area_id   = .data$area_id, 
        parent_area_name = .data$area_name
      )
    areas_model <- areas_model %>% 
      dplyr::left_join(parent_areas, by = "parent_area_id") %>% 
      dplyr::relocate(.data$parent_area_name, .after = .data$parent_area_id)
  }

  # create skeleton dataset with row for every unique area_id, area_name,
  # space, year and circ_age
  out <- tidyr::crossing(
    dplyr::select(areas_model, -dplyr::contains("parent_area")),
    "year"     = seq(start_year, end_year, by = 1),
    "circ_age" = c(0, seq_len(max(survey_circumcision$circ_age, na.rm = TRUE)))
  ) %>%
    # Get time and age variable
    dplyr::mutate(
      time = .data$year - start_year + 1,
      age  = .data$circ_age + 1
    ) %>%
    # Sort dataset
    dplyr::arrange(.data$space, .data$age, .data$time)
  n <- nrow(out)
  
  # Add population data on to merge
  out <- out %>% 
    dplyr::left_join(
      populations %>%
        dplyr::select(
          .data$area_id, 
          .data$year,
          circ_age = .data$age, 
          .data$population
        ),
      by = c("area_id", "year", "circ_age")
    )
  
  # fill in NAs for populations with populations of child areas 
  # TODO: expand to have work for area_level differences > 1
  parent_areas_na_pops <- out %>% 
    filter(is.na(population), area_level < area_lev) %>% 
    distinct(area_id) %>% 
    pull()
  
  if (length(parent_areas_na_pops) != 0) {
    
    message(paste0(
      "Populations missing for ", 
      paste(parent_areas_na_pops, collapse = ", "), 
      ", filling in using populations of child areas"
    ))
    
    # find areas which have parent areas with missing populations
    areas_child_pops <- areas_model %>% 
      dplyr::filter(.data$parent_area_id %in% parent_areas_na_pops) %>% 
      dplyr::select(.data$area_id, .data$parent_area_id, .data$parent_area_name)
    
    # Pull populations for child areas, aggregate to parent area_ids
    missing_pops <- out %>% 
      dplyr::filter(area_id %in% areas_child_pops$area_id) %>% 
      dplyr::left_join(areas_child_pops, by = "area_id") %>% 
      dplyr::mutate(area_id = parent_area_id, area_name = parent_area_name) %>% 
      dplyr::select(-contains("parent")) %>% 
      dplyr::group_by(area_id, year, circ_age, time, age) %>% 
      dplyr::summarise(missing_pops = sum(population), .groups = "drop")
    
    # replace missing populations in skeleton dataset
    out <- out %>% 
      dplyr::left_join(
        missing_pops, 
        by = c("area_id", "year", "circ_age", "time", "age")
      ) %>% 
      mutate(
        population = ifelse(is.na(population), missing_pops, population)
      ) %>% 
      select(-missing_pops)
  }
   
  # give error if there are NAs in populations
  stopifnot(all(!is.na(out$population))) 
  
  # give warning about duplicated rows from left_join
  if (n != nrow(out)) {
    message(paste0(
      "Some rows duplicated by left-joining in populations, function may ",
      "fail below"
    ))
  }

  # Fail if still missing pops, as will lead to more obscure matrix errors
  stopifnot(!all(is.na(out$population)))

  # Add `space` to survey_circumcision observations
  survey_circumcision <- survey_circumcision %>%
    dplyr::left_join(
      dplyr::select(areas_model, .data$area_id, .data$space),
      by = "area_id"
    )
  stopifnot(!is.na(survey_circumcision$space))

  # Obtain N person years
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

  # calculate empirical agetime hazard matrices for different circumcision
  # types, and take column sums (i.e. N empirical circs for each "type):
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

  # add to out:
  out[, empirical_circ_cols] <- agetime_hazard_matrices
  return(out)
}


#' @title Use model shell dataset to estimate empirical circumcision rates
#' @description  Takes the shell dataset with a row for every unique area ID,
#' area name, year and circumcision age in survey data outputed by
#' \link[threemc]{create_shell_dataset} and returns the empirical circumcision
#' rates for each row, aggregated to age groups from single ages. Also converts
#' from wide format to long format.
#' @inheritParams create_shell_dataset
#' @param out Shell dataset outputted by \link[threemc]{create_shell_dataset}
#' @param age_groups Age groups to aggregate by, Default:
#' c("0-4",   "5-9",   "10-14", "15-19",
#'   "20-24", "25-29", "30-34", "35-39",
#'   "40-44", "45-49", "50-54", "54-59",
#'   "15-24", "10-24", "15-29", "10-29",
#'   "15-39", "10-39", "15-49", "10-49"
#' )
#' @seealso
#'  \code{\link[threemc]{create_shell_dataset}}
#' @rdname threemc_empirical_rates
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr %>%
threemc_empirical_rates <- function(out,
                                    areas,
                                    area_lev,
                                    populations,
                                    age_groups = c(
                                      # 5-year age groups from 0 to 60
                                      "0-4", "5-9", "10-14", "15-19",
                                      "20-24", "25-29", "30-34", "35-39",
                                      "40-44", "45-49", "50-54", "54-59",
                                      # other, wider age groups of interest
                                      "0+", "10+", "15+",
                                      "15-24", "10-24", "15-29", "10-29",
                                      "15-39", "10-39", "15-49", "10-49"
                                    )) {

  # wide formatted areas, for changing area levels later
  areas_wide <- areas %>%
    dplyr::select(
      dplyr::contains("area"), .data$space
    ) %>%
    threemc::spread_areas()

  results <- out %>%
    # calculate MC as MC + MMC + TMC
    dplyr::mutate(
      obs_mc = dplyr::case_when(
        !is.na(.data$obs_mmc) &
          !is.na(.data$obs_tmc) ~ .data$obs_mc +
          .data$obs_mmc +
          .data$obs_tmc,
        !is.na(.data$obs_mmc) ~ .data$obs_mc + .data$obs_mmc,
        !is.na(.data$obs_tmc) ~ .data$obs_mc + .data$obs_tmc,
        TRUE ~ .data$obs_mc
      )
    ) %>%
    # pivot empirical person year columns to the one column
    tidyr::pivot_longer(
      cols = .data$obs_mmc:.data$obs_mc,
      names_to = "type",
      values_to = "mean"
    ) %>%
    # Only keep required columns
    dplyr::select(
      .data$area_id:.data$age, .data$type, .data$N, .data$mean
    ) %>%
    # Calculate empirical rates
    dplyr::mutate(
      mean = ifelse(
        .data$mean == 0 | .data$N == 0, 0, .data$mean / .data$N
      )
    ) %>%
    # remove person years column
    dplyr::select(-.data$N) %>%
    # remove 0s, so taking log doesn't give -Inf
    dplyr::mutate(mean = ifelse(.data$mean == 0, NA_real_, .data$mean))

  # only keep relevant columns in populations for left_join
  populations_append <- populations %>%
    dplyr::select(
      names(results)[names(results) %chin% names(populations)],
      .data$population,
      # don't join by area_name, in case character encoding causes errors
      -dplyr::matches("area_name")
    )

  # join in region populations
  results <- dplyr::left_join(results, populations_append)

  if (any(is.na(results$population))) {
    message(paste0(
      "NA populations found, check that population start year is >= ",
      "minimum year in out"
    ))
  }

  # Add parent areas
  results <- combine_areas(
    # for now, ignore results for lower area levels
    dplyr::filter(results, .data$area_level == area_lev),
    areas_wide,
    area_lev,
    add_keep_cols = "mean",
    join = FALSE,
    fill = TRUE
  )

  # ensure there is no duplication
  results <- results %>%
    dplyr::bind_rows() %>%
    dplyr::distinct() %>%
    dplyr::group_split(.data$area_level, .keep = TRUE)

  # aggregate to population-weighted age groups
  results <- aggregate_sample_age_group(
    results,
    aggr_cols  = c("area_id", "area_name", "year", "type"),
    num_cols   = "mean",
    age_groups = age_groups
  )

  # Merge regional information on the dataset
  merge_empirical_rates <- function(.data) {
    merge_area_info(.data, sf::st_drop_geometry(areas)) %>%
      dplyr::mutate(
        iso3 = substr(.data$area_id, 0, 3),
        type = dplyr::case_when(
          # match type with that outputted by threemc_aggregate
          .data$type == "obs_mc" ~ "MC probability",
          .data$type == "obs_mmc" ~ "MMC probability",
          .data$type == "obs_tmc" ~ "TMC probability"
        )
      ) %>%
      dplyr::relocate(.data$iso3) %>%
      dplyr::relocate(dplyr::contains("age"), .after = .data$year)
  }

  return(merge_empirical_rates(results))
}
