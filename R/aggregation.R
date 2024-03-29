#### Main Function ####

#' @title Produce Population Weighted Aggregated Samples for All Area Levels
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
#' @param probs Percentiles to provide quantiles at. Set to NULL to skip
#' computing quantiles.
#' @param ... Further arguments to internal functions.
#' @return \code{data.frame} with samples aggregated by \code{aggr_cols} and
#' weighted by population.
#' @importFrom dplyr %>%
#' @importFrom data.table %chin% .N
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
  probs = c(0.025, 0.5, 0.975),
  ...
  ) {
  
  # produce error if fit does not contain samples
  stopifnot("sample" %in% names(fit))
  
  # produce error if age_var not correctly specified
  stopifnot(length(age_var) == 1 && age_var %in% c("age", "age_group"))
  
  # copy data.table datasets so destructive operations don't change global vars
  .data <- data.table::copy(.data)
  areas <- data.table::copy(areas)
  
  # global bindings for data.table non-standard evaluation
  space <- age <- NULL
  
  #### Prepare location/shapefile information ####
  if (inherits(areas, "sf")) {
    areas <- sf::st_drop_geometry(areas)
  }
  
  # fill missing pops for historical/future years not in data, if required
  start_year <- min(.data$year)
  end_year <- max(.data$year)
  min_pop_year <- min(populations$year)
  max_pop_year <- max(populations$year)
  if (start_year < min_pop_year || end_year > max_pop_year) {
    populations <- fill_downup_populations(
      populations, 
      start_year, 
      end_year,
      min_pop_year,
      max_pop_year
    )
  }
  
  # Add a unique identifier within Admin code and merging to boundaries
  data.table::setDT(areas)[,  space := seq_len(.N), by = "area_level"]

  # wide formatted areas, for changing area levels later
  areas_wide <- areas[,
      c("area_id", "area_name", "parent_area_id", "area_level", "space")
  ]
  areas_wide <- spread_areas(areas_wide, space = FALSE)

  # Model with Probability of MC
  .data$model <- "No program data"
  fit_no_prog <- fit # need to load model with programme data as well
  
  # Take sample matrix rows that are kept in .data from original skeleton data
  if ("n" %in% names(.data)) {
    fit_no_prog$sample <- lapply(fit_no_prog$sample, function(x) x[.data$n, ])
    .data$n <- NULL
  } 
  # throw error if number of rows in results does not equal sample number
  stopifnot(nrow(.data) == nrow(fit_no_prog$sample$haz))

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
  num_cols <- c(paste0("samp_", seq_len(N)))
  if (age_var == "age") {
    .data <- aggregate_sample(.data, num_cols = num_cols, ...)
  } else {
    .data <- aggregate_sample_age_group(.data, num_cols = num_cols, ...)
  }

  # additional aggregations to perform for prevalence
  if (type == "prevalence" &&
      !is.null(prev_year) &&
      prev_year %in% .data$year) {

    # calculate change in prevalence since prev_year
    data_change_prev_year <- prevalence_change(
      .data,
      spec_year = prev_year
    )

    # Get number of people circumcised
    data_n <- n_circumcised(.data)

    # add "change from prev_year" and "n_circumcised" .data to .data
    .data <- data.table::rbindlist(
      list(.data, data_change_prev_year, data_n),
      use.names = TRUE
    )
  }
  .data <- posterior_summary_fun(.data, probs = probs)

  # Merge regional information on the dataset (i.e. parent area info) & return
  .data <- merge_area_info(.data, areas)
}

#### prepare_sample_data ####

#' @title Pull N Circumcision Samples From TMB Fit
#' @description Function to pull samples for various summaries inferred about
#' circumcision in a given area (prevalence, probability or incidence,
#' respectively).
#' @inheritParams threemc_aggregate
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

  # global bindings for data.table non-standard evaluation
  model <- NULL

  #### Load rates from survival model ####

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
      # fill needed for no type models as tmpx_i may have all NAs
      data.table::rbindlist(use.names = TRUE, fill = TRUE) %>%
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
    tmp <- data.table::merge.data.table(
      tmp, populations_append, all.x = TRUE
    ) %>%
      dplyr::relocate(.data$population, .before = .data$samp_1)

    # filter out NA pops, return appropriate message
    if (any(is.na(tmp$population)) == TRUE) {
      # original number of rows and years
      n1 <- nrow(tmp)
      years <- unique(tmp$year)
      tmp <- dplyr::filter(tmp, !is.na(.data$population))
      n2 <- nrow(tmp)
      if (n2 == 0) stop("No populations present in data")
      years_pop <- unique(tmp$year)
      years_no_pop <- years[!years %in% years_pop]
      message(paste0("Missing population for ", n1 - n2, " records\n"))
      message(paste0("Years removed: ", paste(years_no_pop, collapse = ", ")))
    }

    return(tmp)
  }

  # Model with Probability of MC with no program data (only surveys)
  if (!is.null(no_prog_results)) {
    tmp1 <- data.table::copy(
      data.table::setDT(no_prog_results)[, model := "No program data"]
    )
    tmp1 <- append_fun(
      as.data.frame(tmp1), no_prog_tmb_fit, populations, type = type
    )
  } else {
    tmp1 <- NULL
  }

  # Model with Probability of MC with both programme and survey data
  if (!is.null(prog_results)) {
    tmp2 <- data.table::copy(
      data.table::setDT(prog_results)[, model := "With program data"]
    )
    tmp2 <- append_fun(
      as.data.frame(tmp2), prog_tmb_fit, populations, type = type
    )
  } else {
    tmp2 <- NULL
  }

  # Append together
  return(data.table::rbindlist(list(tmp1, tmp2)))
}

#### aggregate sample ####

#' @title Produce Population Weighted Aggregated Samples
#' @description Aggregate by area, year, age and type (weighted by population),
#' and convert to a percentage/probability.
#' @param .data \code{data.frame} including area populations, with
#' un-aggregated samples.
#' @param aggr_cols Columns to aggregate samples by, Default:
#' c("area_id", "area_name", "year", "age", "age_group", "model", "type")
#' @param num_cols `numeric` columns to aggregate.
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
                             ),
                             num_cols,
                             ages = 0:60, 
                             ...) {

  # global bindings for data.table non-standard evaluation
  .SD <- population <- age <- NULL

  # convert .data from list to data.frame, if required
  if (!inherits(.data, "data.frame")) {
    .data <- data.table::rbindlist(
      .data,
      use.names = TRUE,
      fill = TRUE,
      ...
    )
  }
  
  # only keep specified ages
  .data <- .data[age %in% ages, ]
  
  # remove population from num_cols, if present 
  num_cols <- num_cols[!num_cols == "population"]

  # Multiply by population to population weight
  .data <- data.table::setDT(.data)[,
    (num_cols) := lapply(.SD, function(x) x * population),
    .SDcols = num_cols
  ]

  # ensure aggregation columns are in the data
  aggr_cols <- aggr_cols[aggr_cols %chin% names(.data)]

  # summarise num_cols (along with `population`) by aggr_cols
  sd_cols <- grep(
    paste(c("population", num_cols), collapse = "|"), names(.data)
  )
  .data <- .data[,
    lapply(.SD, sum, na.rm = TRUE),
    by = c(aggr_cols),
    .SDcols = sd_cols
  ]

  # divide by population to population weight
  return(.data[,
          (num_cols) := lapply(.SD, function(x) {
            # data.table::fifelse(
            #   grepl("performed", type), x, x / population
            # )
            x / population
          }),
          .SDcols = num_cols
  ])
}

#### aggregate_sample_age_group ####

#' @title Produce Population Weighted Aggregations for Age Groups
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
      "0-14",  "10-24", "15-24", "10-29", "15-29",
      "10-39", "15-39", "10-49", "15-49", "30-49"
    )
) {

  # global bindings for `data.table` non-standard evaluation
  .SD <- population <- type <- age_group <- NULL

  # bind list objects if required
  results <- results_list
  if (!inherits(results, "data.frame")) {
    results <- data.table::rbindlist(results, use.names = TRUE)
  }

  # populations, aggr_cols and num_cols should all be present in results
  stopifnot(all(c("population", aggr_cols, num_cols) %chin% names(results)))

  # convert to data.table
  data.table::setDT(results)
  # avoid duplication
  results <- unique(results)
  # Multiply num_cols by population to population weight
  results[,
          (num_cols) := lapply(.SD, function(x) x * population),
          .SDcols = num_cols
        ]
  # have num_cols at the end
  other_names <- names(results)[!names(results) %chin% num_cols]
  results <- data.table::setcolorder(results, c(other_names, num_cols))

  # create data frame matching age groups to ages within
  age_group_df <- data.table::rbindlist(
    lapply(age_groups, match_age_group_to_ages, max_age = max(results$age))
  )

  # left join in  `age_group` col to results, remove `age` col
  results <- data.table::merge.data.table(
    results,
    age_group_df,
    by = "age",
    all.x = TRUE,
    # allow duplication of records falling in multiple age groups
    allow.cartesian = TRUE
  )[, !c("age")]
  
  # remove missing age groups
  results <- results[!is.na(age_group), ]

  if (!"age_group" %chin% aggr_cols) aggr_cols <- c(aggr_cols, "age_group")

  # aggregate sample for unique combination of `aggr_cols`
  sd_cols <- c("population", num_cols)
  results <- results[,
                     lapply(.SD, sum, na.rm = TRUE),
                     by = c(aggr_cols),
                     .SDcols = sd_cols
             ]

  # create dummy type column, if required
  if (!"type" %chin% names(results)) results$type <- "dummy"

  # Divide by population to population weight
  results[,
    (num_cols) := lapply(.SD, function(x) x / population), .SDcols = num_cols
  ]

  # remove dummy column
  if (all(results$type == "dummy")) results <- results[, -c("type")]

  # return results
  return(results)
}


#### prevalence_change ####

#' @title Calculate Change in Prevalence/Coverage From a Given Year
#' @description Function to calculate change in prevalence/coverage from a
#' given year for all other years.
#' @param results results for several years.
#' @param spec_year Year to calculate change in prevalence/coverage from within
#' \code{results}.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname prevalence_change
#' @keywords internal
prevalence_change <- function(results, spec_year) {

  # global bindings for data.table non-standard evaluation
  year <- value <- prev_value <- type <- NULL

  results <- data.table::setDT(results)

  # pull samples from coverage in chosen year
  spec_year_results <- results[
    year == spec_year,
    !c("year", "population")
  ]
  spec_year_results <- data.table::melt(
    spec_year_results,
    measure = patterns("^samp"),
    value.name = "prev_value"
  )

  if (nrow(spec_year_results) == 0) {
    stop("Please choose a different year to draw comparisons with")
  }

  # filter for years above comparison year and pivot longer
  results <- data.table::melt(
    results[year > spec_year], measure = patterns("^samp")
  )

  # join into spec_year_results for corresponding categorical vars
  results <- data.table::merge.data.table(
    results, spec_year_results, all.x = TRUE
  )
  # subtract value from comparison year from each year
  results <- results[,
                     value := value - prev_value
  ][,
    !c("prev_value")
  ]

  # pivot back to wide format
  results <-  data.table::dcast(results, ... ~ variable, value.var = "value")[,
                type := paste0("Change in ", type, " from ", spec_year)
              ]
}

#### n_circumcised ####

#' @title Calculate Number of People Circumcised
#' @description Calculate number of people circumcised (as well as unmet need).
#' @param results Results with samples for number of circumcisions performed
#' in each region.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname n_circumcised
#' @keywords internal
n_circumcised <- function(results) {

  # global bindings for `data.table` non-standard evaluation
  type <- population <- NULL

  # Get number of circumcised men
  results <- data.table::copy(data.table::setDT(results))[,
                               type := paste0(
                                 "Number circumcised (",
                                 stringr::str_remove(type, " coverage"),
                                 ")"
                               )
                              ]

  # split by type (uses split.data.table method)
  n_circ <- split(results, by = "type")

  # get circumcised population by type
  sd_cols <- grep("samp_", names(results))
  n_circ_type <- lapply(n_circ, function(x) {
    x <- data.table::setDT(x)[,
                         (sd_cols) := lapply(.SD, function(y) y * population),
                         .SDcols = sd_cols
                         ]
    return(x)
  })

  # also calculate unmet need for total circumcisions
  unmet_need <- data.table::copy(results)[
    type == "Number circumcised (MC)"
  ][,
    (sd_cols) := lapply(.SD, function(x) population * (1 - x)),
    .SDcols = sd_cols
  ][,
    type := "Unmet need"
    ]

  # ensure unmet need does not go below 0
  samp_cols <- names(unmet_need)
  samp_cols <- samp_cols[grepl("samp_", samp_cols)]
  unmet_need[,
    (samp_cols) := lapply(.SD, function(x) data.table::fifelse(x < 0, 0, x)),
    .SDcols = samp_cols
  ]

  n_circ_type[[length(n_circ_type) + 1]] <- unmet_need

  # Append together
  (data.table::rbindlist(n_circ_type, use.names = TRUE))
}

#### posterior_summary_fun ####

#' @title Calculate Summary Statistics From Samples
#' @description Takes samples and calculates summary statistics (mean, standard
#' deviation, and quantiles (if desired)).
#' @inheritParams threemc_aggregate
#' @param .data \code{data.frame} with samples to be summarised.
#' @importFrom dplyr %>%
#' @importFrom rlang :=
#' @rdname posterior_summary_fun
#' @keywords internal
posterior_summary_fun <- function(.data, probs = c(0.025, 0.5, 0.975)) {

  # global bindings for data.table non-standard evaluation
  . <- value <- NULL

  # ensure probabilities for quantiles are ordered
  probs <- sort(probs)

  .data <- data.table::setDT(.data)

  # ensure numeric columns (apart from age) are after categorical
  num_cols <- names(.data)
  num_cols <- c("population", num_cols[grep("samp_", num_cols)])
  data.table::setcolorder(.data)

  # pull locations of columns to "group by"
  id_cols <- seq_along(names(.data)[!grepl("samp", names(.data))])

  # use data.table as this can be quite slow for larger countries
  .data <- data.table::setDT(.data)

  # pivot to calculate row means and sds of samples for each stratification
  samp_cols <- grep("samp", names(.data))
  .data_long <- data.table::melt(.data,
    id.vars = id_cols,
    measure.vars = samp_cols
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

    return(data.table::merge.data.table(.data, quantiles))
  } else {
    return(.data)
  }
}


#### merge_area_info ####

#' @title Merge Regional Information On Dataset
#' @description Merge regional information on the dataset
#' (i.e. parent area info).
#' @param results \code{data.frame} you wish to merge shapefiles with.
#' @inheritParams threemc_aggregate 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname merge_area_info
#' @keywords internal
merge_area_info <- function(results, areas) {

  # global bindings for `data.table` non-standard evaluation
  . <- area_id <- area_name <- area_level <- ..rm_cols <- ..keep_cols <- NULL

  if (!inherits(results, "data.table")) results <- data.table::setDT(results)
  if (!inherits(areas, "data.table")) areas <- data.table::setDT(areas)

  # remove area_name & area_level to avoid mishaps when joining
  rm_cols <- c("area_name", "area_level")
  rm_cols <- rm_cols[rm_cols %chin% names(results)] # avoids warning
  results <- results[, !..rm_cols]

  # Add region information
  results <- data.table::merge.data.table(
    x     = results,
    y     = areas[, .(area_id, area_name, area_level)],
    by    = "area_id", 
    all.x = TRUE
  )

  # Merge regional information on the dataset (i.e. parent area info)
  # use parent area names from areas, if one of parent_names is in results
  parent_names <- c("parent_area_id", "parent_area_name")
  keep_cols <- c("area_id", parent_names)
  keep_cols <- keep_cols[keep_cols %in% names(areas)]
  areas_join <- areas[, ..keep_cols]
  # add parent_area_name using area_name from areas, if missing
  if (!"parent_area_name" %in% names(areas_join)) {
    areas_join <- data.table::merge.data.table(
      x = areas_join, 
      y = data.table::setnames(
        areas[, .(area_id, area_name)],
        old = c("area_id", "area_name"),
        new = parent_names
      ),
      all.x = TRUE
    )
  }
  results <- data.table::merge.data.table(
    x     = results,
    y     = areas_join,
    all.x = TRUE
  )   

  # relocate area cols to the start of the dataframe
  area_cols <- c(
    "area_name", "area_id", "area_level", "parent_area_id", "parent_area_name"
  )
  area_cols <- area_cols[area_cols %in% names(results)]
  data.table::setcolorder(
    results,
    c(area_cols, dplyr::setdiff(names(results), area_cols))
  )
}
