#### Main Function ####

#' @title Produce Population weighted Aggregated Samples
#' @description Aggregate by area, year, age and type (weighted by population),
#' and convert to a percentage/probability.
#' @param .data \code{data.frame} of unaggregated modelling results.
#' @param fit \code{TMB} list containing model parameters, nested list of
#' `samples` for the (cumulative) incidence and hazard rate of circumcision for
#' the region(s) in question.
#' @param areas `sf` shapefiles for specific country/region.
#' @param populations \code{data.frame} containing populations for each
#' region in tmb fits.
#' @param age_var Determines whether you wish to aggregate by discrete ages or
#' age groups (0-4, 5-9, 10-14, and so on).
#' @param type Determines which aspect of MC in the regions in question you wish
#' to aggregate for. Can be one of "probability", "incidence" or "prevalence".
#' @param area_lev PSNU area level for specific country.
#' @param N Number of samples to be generated, Default: 100
#' @param prev_year If type == "prevalence", choose year to compare prevalence
#' with.
#' @param ... Further arguments to internal functions.
#' @return \code{data.frame} with samples aggregated by \code{aggr_cols} and
#' weighted by population.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname threemc_aggregate
#' @export
threemc_aggregate <- function(
  # datasets
  .data, fit, areas, populations, # datasets
  # options
  age_var = c("age", "age_group"),
  type = c("probability", "incidence", "prevalence"),
  area_lev,
  N = 100,
  prev_year = 2008,
  ...
  ) {

  #### Preparing location/shapefile information ####
  if (inherits(areas, "sf")) {
    areas <- sf::st_drop_geometry(areas)
  }

  areas <- areas %>%
    # Add a unique identifier within Admin code and merging to boundaries
    dplyr::group_by(.data$area_level) %>%
    dplyr::mutate(space = dplyr::row_number()) %>%
    dplyr::ungroup()

  # wide formatted areas, for changing area levels later
  areas_wide <- areas %>%
    dplyr::select(dplyr::all_of(
        c("area_id", "area_name", "parent_area_id", "area_level", "space")
    )) %>%
    spread_areas(space = FALSE)

  # Model with Probability of MC
  .data$model <- "No program data"
  fit_no_prog <- fit # need to load model with programme data as well

  #### Load rates from survival model ####
  .data <- prepare_sample_data(
    N = N,
    populations = populations,
    no_prog_results = .data,
    no_prog_tmb_fit = fit_no_prog,
    type = type
  )

  #### Prepare .data for output ####

  # collect .data for lower area hierarchies by joining higher area
  # hierarchies (do so until you reach "area 0")
  .data <- combine_areas(.data, areas_wide, area_lev, join = FALSE)

  # aggregate samples for each individual age or age group
  if (age_var == "age") {
    .data <- aggregate_sample(.data, ...)
  } else {
    .data <- aggregate_sample_age_group(
      .data,
      num_cols = c(paste0("samp_", seq_len(N))),
      ...
    )
  }

  # additional aggregations to perform for prevalence
  if (type == "prevalence" && prev_year %chin% .data$year) {
    # calculate change in prevalence since prev_year
    data_change_prev_year <- prevalence_change(
      .data,
      spec_year = prev_year
    )

    # Getting number of people circumcised
    data_n <- n_circumcised(.data)

    # add "change from prev_year" and "n_circumcised" .data to .data
    .data <- rbind(.data, data_change_prev_year, data_n)
  }
  # calculate summary statistics (mean, sd, quantiles)
  .data <- posterior_summary_fun(.data, ...)

  # Merge regional information on the dataset (i.e. parent area info) & return
  .data <- merge_area_info(.data, areas)
}

#### prepare_sample_data ####

#' @title Pull N circumcision samples from TMB fit
#' @description Function to pull samples for various summaries inferred about
#' circumcision in a given area (prevalence, probability or incidence,
#' respectively).
#' @param N Number of samples to be generated, Default: 100
#' @param populations \code{data.frame} containing populations for each
#' region in tmb fits.
#' @param no_prog_results Quantiles of different summaries for different types
#' (Medical, traditional or total) circumcision, for survey data only,
#' Default: NULL
#' @param prog_results Quantiles of different summaries for different types
#' (Medical, traditional or total) circumcision, for survey data and VMMC
#' programme data, Default: NULL
#' @param no_prog_tmb_fit TMB model for Male Circumcision (MC), for survey data
#' only.
#' @param prog_tmb_fit TMB model for Male Circumcision (MC), for survey data and
#' VMMC programme data.
#' @param type Determines which aspect of MC in the regions in question you wish
#' to sample for. Can be one of "probability", "incidence" or "prevalence".
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname prepare_sample_data
#' @keywords internal
prepare_sample_data <- function(N = 100,
                                populations,
                                no_prog_results = NULL,
                                prog_results = NULL,
                                no_prog_tmb_fit,
                                prog_tmb_fit,
                                type) {
  if (is.null(no_prog_results) && is.null(prog_results)) {
    stop("cannot have prog_results == no_prog_results == NULL")
  }
  if (!type %chin% c("probability", "incidence", "prevalence") ||
    length(type) > 1) {
    stop("Please choose a valid type
         (one of 'probability', 'incidence', 'prevalence'")
  }

  if (nrow(populations) == 0) stop("No populations present in data")

  #########################################
  ### Loading rates from survival model ###
  #########################################

  # append samples from different circumcision models together (and pops)
  append_fun <- function(tmp, fit, populations, type) {
    
    # different objects to pull samples from fit, based on desired aggregation
    if (type == "probability") {
      mmc <- "haz_mmc" # medical circumcision
      tmc <- "haz_tmc" # traditional circumcision
      mc <- ifelse("haz" %chin% names(fit$sample), "haz", "haz_mc") # all circ
    } else if (type == "incidence") {
      mmc <- "inc_mmc"
      tmc <- "inc_tmc"
      mmct <- "inc_mmct"
      mc <- ifelse("inc" %chin% names(fit$sample), "inc", "inc_mc")
    } else if (type == "prevalence") {
      mmc <- "cum_inc_mmc"
      tmc <- "cum_inc_tmc"
      mmct <- "cum_inc_mmct"
      mc <- ifelse("cum_inc" %chin% names(fit$sample), "cum_inc", "cum_inc_mc")
    }
    
    # word to be pasted onto the end of circ type below
    if (type == "prevalence") {
      category <- "coverage"
    } else {
      category <- type
    }
    
    # initialise dataframes to store samples for different circ types
    tmpx_1 <- tmpx_2 <- tmpx_3 <- tmpx_4 <- tmpx_5 <- tmpx_6 <- tmp
    # type == "incidence" has 12 different "types"
    if (type == "incidence") {
      tmpx_7 <- tmpx_8 <- tmpx_9 <- tmpx_10 <- tmpx_11 <- tmpx_12 <- tmp
    } else {
      tmpx_7 <- tmpx_8 <- tmpx_9 <- tmpx_10 <- tmpx_11 <- tmpx_12 <- NULL
    }
    
    # Pull samples from data
    # Models with no VMMC data cannot distinguish between MMC-nT and MMC-T
    tmpx_1[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][, 1:N]
    tmpx_1$type <- paste("MMC-nT", category)
    if (tmp$model[1] == "No program data") {
      tmpx_2[, paste0("samp_", 1:N)] <- 0
      tmpx_3[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N]
      tmpx_4[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][, 1:N]
      tmpx_5[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N]
      tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[mc]][, 1:N]
    } else if (tmp$model[1] == "With program data") {
      if (type == "probability") {
        tmpx_2[, paste0("samp_", 1:N)] <- fit$sample$probs[, 1:N] *
          fit$sample[[tmc]][, 1:N]
        tmpx_3[, paste0("samp_", 1:N)] <- (1 - fit$sample$probs[, 1:N]) *
          fit$sample[[tmc]][, 1:N]
        tmpx_4[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][, 1:N] +
          fit$sample$probs[, 1:N] * fit$sample[[tmc]][, 1:N]
        tmpx_5[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N]
        tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[mc]][, 1:N]
      } else {
        tmpx_2[, paste0("samp_", 1:N)] <- fit$sample[[mmct]][, 1:N]
        tmpx_3[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N]
        tmpx_4[, paste0("samp_", 1:N)] <- fit$sample[[mmc]][, 1:N] +
          fit$sample[[mmct]][, 1:N]
        tmpx_5[, paste0("samp_", 1:N)] <- fit$sample[[tmc]][, 1:N] +
          fit$sample[[mmct]][, 1:N]
        tmpx_6[, paste0("samp_", 1:N)] <- fit$sample[[mc]][, 1:N]
      }
    }
    # give appropriate labels to each df
    tmpx_2$type <- paste("MMC-T", category)
    tmpx_3$type <- paste("TMC", category)
    tmpx_4$type <- paste("MMC", category)
    tmpx_5$type <- paste("TMIC", category)
    tmpx_6$type <- paste("MC", category)
    
    # Samples for the number of MCs performed (for incidence)
    if (type == "incidence") {
      tmpx_7 <- tmpx_1
      tmpx_7$type <- "MMC-nTs performed"
      tmpx_8 <- tmpx_2
      tmpx_8$type <- "MMC-Ts performed"
      tmpx_9 <- tmpx_3
      tmpx_9$type <- "TMCs performed"
      tmpx_10 <- tmpx_4
      tmpx_10$type <- "MMCs performed"
      tmpx_11 <- tmpx_5
      tmpx_11$type <- "TMICs performed"
      tmpx_12 <- tmpx_6
      tmpx_12$type <- "MCs performed"
    }
    
    # Append together
    tmp <- as.list(mget(paste0("tmpx_", 1:12))) %>%
      data.table::rbindlist() %>%
      # only keep relevant columns
      dplyr::select(
        .data$area_id, .data$area_name,
        .data$year, .data$age,
        .data$type, .data$model,
        dplyr::contains("samp_")
      )
    
    # only keep relevant columns in populations for left_join
    populations_append <- populations %>%
      dplyr::select(
        dplyr::all_of(names(tmp)[names(tmp) %chin% names(populations)]),
        .data$population,
        # don't join by area_name, in case character encoding etc causes errors
        -dplyr::matches("area_name")
      )
    
    # join with populations
    tmp <- data.table::merge.data.table(tmp, populations_append) %>% 
      dplyr::relocate(.data$population, .before = .data$samp_1)
    
    # filter out na populations, with an appropriate message
    if (any(is.na(tmp$population)) == TRUE) {
      n1 <- nrow(tmp)
      tmp <- dplyr::filter(tmp, !is.na(.data$population))
      n2 <- nrow(tmp)
      if (n2 == 0) stop("No populations present in data")
      message(paste0("Missing population for ", n1 - n2, " records"))
    }
    
    return(tmp)
  }
  
  # Model with Probability of MC with no program data (only surveys)
  if (!is.null(no_prog_results)) {
    tmp1 <- dplyr::mutate(no_prog_results, model = "No program data")
    tmp1 <- append_fun(tmp1, no_prog_tmb_fit, populations, type = type)
  } else {
    tmp1 <- NULL
  }

  # Model with Probability of MC with both programme and survey data
  if (!is.null(prog_results)) {
    tmp2 <- dplyr::mutate(prog_results, model = "With program data")
    tmp2 <- append_fun(tmp2, prog_tmb_fit, populations, type = type)
  } else {
    tmp2 <- NULL
  }

  # Append together
  return(data.table::rbindlist(list(tmp1, tmp2)))
}

#### aggregate sample ####

#' @title Produce Population weighted Aggregated Samples
#' @description Aggregate by area, year, age and type (weighted by population),
#' and convert to a percentage/probability.
#' @param .data \code{data.frame} including area populations, with
#' un-aggregated samples.
#' @param aggr_cols Columns to aggregate samples by, Default:
#' c("area_id", "area_name", "year", "age", "age_group", "model", "type")
#' @param ... Further arguments passed to \code{data.table::rbindlist}.
#' @return \code{data.frame} with samples aggregated by \code{aggr_cols} and
#' weighted by population.
#' @importFrom dplyr %>%
#' @importFrom data.table .SD
#' @rdname aggregate_sample
#' @keywords internal
aggregate_sample <- function(.data,
                             aggr_cols = c(
                               "area_id", "area_name", "year",
                               "age", "age_group", "model", "type"
                             ), ...) {

  # global bindings for data.table non-standard evaluation
  .SD <- NULL

  # convert .data from list to data.frame, if required
  if (!inherits(.data, "data.frame")) {
    .data <- as.data.frame(data.table::rbindlist(
      .data,
      use.names = TRUE,
      fill = TRUE,
      ...
    ))
  }

  # ensure aggregation columns are in the data
  aggr_cols <- aggr_cols[aggr_cols %chin% names(.data)]

  # Multiplying by population to population weight
  .data <- .data %>%
    dplyr::mutate(dplyr::across(dplyr::contains("samp_"), ~ . * population))

  # summarise samples by aggr_cols
  .data <- data.table::setDT(.data)[,
    lapply(.SD, sum, na.rm = TRUE),
    by = c(aggr_cols),
    .SDcols = c("population", paste0("samp_", c(1:100)))
  ]
  # divide by population to population weight
  .data <- as.data.frame(.data) %>%
    dplyr::mutate(dplyr::across(dplyr::contains("samp_"), ~ . / population))

  return(.data)
}

#### aggregate_sample_age_group ####

#' @title Produce Population weighted Aggregations for Age Groups
#' @description Aggregate specified `numeric` columns by population-weighted
#' age groups (rather than single year ages), split by specified categories.
#' @param results_list list of \code{data.frame}s outputted by
#' \code{\link[threemc]{combine_areas}} with \code{join = FALSE}, including
#' area populations, with un-aggregated samples.
#' @param aggr_cols Columns to aggregate samples by, Default:
#' c("area_id", "area_name", "year", "model", "type")
#' @param num_cols `numeric` columns to aggregate.
#' @param age_groups Age groups to aggregate by, Default:
#' c("0-4",   "5-9",   "10-14", "15-19", "20-24", "25-29",
#' "30-34", "35-39", "40-44", "45-49", "50-54", "54-59",
#' "0+",    "10+",   "15+",   "15-24", "10-24", 15-29",
#' "10-29", "15-39", "10-39", "15-49", "10-49")
#' @return \code{data.frame} with samples aggregated by \code{aggr_cols} and
#' weighted by population.
#' @seealso
#'  \code{\link[threemc]{combine_areas}}
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname aggregate_sample_age_group
#' @keywords internal
aggregate_sample_age_group <- function(
    results_list,
    aggr_cols  = c("area_id", "area_name", "year", "model", "type"),
    num_cols,
    age_groups = c(
      # five-year age groups
      "0-4",   "5-9",   "10-14", "15-19", "20-24", "25-29",
      "30-34", "35-39", "40-44", "45-49", "50-54", "54-59",
      # age groups with only minimum cut-off
      "0+", "10+", "15+",
      # other, wider age groups of interest
      "10-24", "15-24", "10-29", "15-29",
      "10-39", "15-39", "10-49", "15-49"
    )
) {
  
  # global bindings for `data.table` non-standard evaluation
  .SD <- NULL
  
  results <- data.table::rbindlist(results_list) %>% 
    # avoid duplication
    dplyr::distinct() %>% 
    # Multiply by population to population weight
    dplyr::mutate(
      dplyr::across(dplyr::any_of(num_cols), ~ . * .data$population)
    ) %>%
    dplyr::relocate(dplyr::any_of(num_cols), .after = dplyr::everything())
  
  # create data frame matching age groups to ages within 
  age_group_df <- data.table::rbindlist(
    lapply(age_groups, match_age_group_to_ages, max_age = max(results$age))
  )
  
  # left join in  `age_group` col to results, remove `age` col
  results <- dplyr::left_join(results, age_group_df, by = "age") %>% 
    dplyr::select(-.data$age)
  
  if (!"age_group" %chin% aggr_cols) aggr_cols <- c(aggr_cols, "age_group")
  
  # aggregate sample for unique combination of `aggr_cols`
  results <- data.table::setDT(results)[,
                                        lapply(.SD, sum, na.rm = TRUE),
                                        by = c(aggr_cols),
                                        .SDcols = c("population", num_cols)
  ]
  
  # Multiply by population to population weight
  # (don"t do this for "N performed", if present)
  results <- results %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::any_of(num_cols), ~ ifelse(
          grepl("performed", type), ., . / population
        )
      )
    )
  
  return(results)
}


#### prevalence_change ####

# function to get change in prevalence/coverage from a given year
#' @title Calculate Change in Prevalence/Coverage from a given year.
#' @description Function to calculate change in prevalence/coverage from a
#' given year for all other years.
#' @param results \code{results} results for several years.
#' @param spec_year Year to calculate change in prevalence/coverage from within
#' \code{results}.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname prevalence_change
#' @keywords internal
prevalence_change <- function(results, spec_year) {

  # pull samples from coverage in chosen year
  spec_year_results <- results %>%
    dplyr::filter(.data$year == spec_year) %>%
    dplyr::select(-c(.data$year, .data$population)) %>%
    tidyr::pivot_longer(dplyr::contains("samp_"), values_to = "prev_value")

  if (nrow(spec_year_results) == 0) {
    stop("Please choose a different year to draw comparisons with")
  }

  # join into spec_year_results for corresponding categorical vars & subtract
  results %>%
    dplyr::filter(.data$year > spec_year) %>%
    tidyr::pivot_longer(dplyr::contains("samp_")) %>%
    dplyr::left_join(spec_year_results) %>%
    dplyr::mutate(value = .data$value - .data$prev_value) %>%
    dplyr::select(-.data$prev_value) %>%
    tidyr::pivot_wider(names_from = "name", values_from = "value") %>%
    dplyr::mutate(type = paste0("Change in ", .data$type, " from ", spec_year))
}

#### n_circumcised ####

#' @title Calculate number of people circumcised
#' @description Calculate number of people circumcised (as well as unmet need).
#' @param results Results with samples for number of circumcisions performed
#' in each region.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname n_circumcised
#' @keywords internal
n_circumcised <- function(results) {

  # Getting number of circumcised men
  n_circ <- dplyr::group_split(results, .data$type)

  # get circumcised population by type
  n_circ_type <- lapply(n_circ, function(x) {
    x %>%
      dplyr::mutate(
        dplyr::across(dplyr::contains("samp_"), ~ . * population),
        type = paste0(
          "Number circumcised (",
          stringr::str_remove(.data$type, " coverage"),
          ")"
        )
      )
  })
  # also calculate unmet need
  n_circ_type[[length(n_circ_type) + 1]] <- results %>%
    dplyr::filter(.data$type == "MC coverage") %>%
    dplyr::mutate(
      dplyr::across(dplyr::contains("samp_"), ~ population * (1 - .)),
      type = "Unmet need"
    )

  # Append together
  (as.data.frame(data.table::rbindlist(n_circ_type, use.names = TRUE)))
}

#### posterior_summary_fun ####

# Getting summary statistics
#' @title Calculate summary statistics from Samples
#' @description Takes samples and calculates summary statistics (mean, standard
#' deviation, and quantiles (if desired)).
#' @param .data \code{data.frame} with samples to be summarised.
#' @param probs Percentiles to provide quantiles at. Set to NULL to skip
#' computing quantiles.
#' @importFrom dplyr %>%
#' @importFrom rlang :=
#' @rdname posterior_summary_fun
#' @keywords internal
posterior_summary_fun <- function(.data, probs = c(0.025, 0.5, 0.975)) {

  # global bindings for data.table non-standard evaluation
  . <- value <- NULL

  probs <- sort(probs)

  # ensure numeric columns are after categorical
  .data <- .data %>%
    dplyr::relocate(
      .data$population | dplyr::contains("samp_"),
      .after = dplyr::everything()
    )

  # pull locations of columns to "group by"
  id_cols <- seq_along(names(.data)[!grepl("samp", names(.data))])

  # use data.table as this can be quite slow for larger countries
  if (!inherits(.data, "data.table")) .data <- data.table::setDT(.data)

  # pivot to calculate row means and sds of samples for each stratification
  .data_long <- data.table::melt(.data,
    id.vars = id_cols,
    measure.vars = c(paste0("samp_", 1:100))
  )
  .data <- .data_long[,
    "."
    (mean = mean(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE)),
    keyby = c(names(.data)[id_cols])
  ] # group by all categories]

  # calculate median and CI
  if (!is.null(probs)) {
    quantiles <- .data_long[, {
        quantiles <- stats::quantile(value,
          probs,
          na.rm = TRUE,
          names = FALSE
        )
        ":="
        list(
          lower  = quantiles[1],
          median = quantiles[2],
          upper  = quantiles[3]
        )
      },
      keyby = c(names(.data)[id_cols]),
    ]

    return(as.data.frame(merge(.data, quantiles)))
  } else {
    return(.data)
  }
}


#### merge_area_info ####

#' @title Merge Regional Informatoin on Dataset
#' @description Merge regional information on the dataset
#' (i.e. parent area info).
#' @param results \code{data.frame} you wish to merge shapefiles with.
#' @param areas \code{sf} shapefiles for specific country/region.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname merge_area_info
#' @keywords internal
merge_area_info <- function(results, areas) {

  # Merging regional information on the dataset (i.e. parent area info)
  results <- results %>%
    # remove area_name & area_level to avoid mishaps when joining
    dplyr::select(
      -c(dplyr::matches("area_name"), dplyr::matches("area_level"))
    ) %>%
    # Adding region information
    dplyr::left_join(
      (areas %>%
        dplyr::select(.data$area_id:.data$area_level)),
      by = c("area_id")
    ) %>%
    dplyr::relocate(.data$area_level, .after = "area_name") %>%
    dplyr::left_join(
      (areas %>%
        dplyr::select(
          parent_area_id = .data$area_id,
          parent_area_name = .data$area_name
        )),
      by = "parent_area_id"
    ) %>%
    dplyr::relocate(dplyr::contains("area"))
}
