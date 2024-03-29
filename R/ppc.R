#' @title Posterior Predictive Distribution and checks on OOS survey
#' @description Aggregate specified `numeric` columns by population-weighted
#' age groups (rather than single year ages), split by specified categories.
#' @param fit Fit object returned by \code{naomi::sample_tmb}, which includes,
#' among other things, the optimised parameters and subsequent sample for our
#' TMB model.
#' @param out Results of model fitting (at specified model `area_lev`) 
#' outputted by \link[threemc]{compute_quantiles}.
#' @param survey_circumcision_test `survey_circumcision` dataset loaded with 
#' \link[threemc]{read_circ_data}. *Do not preprocess with*
#' *\link[threemc]{prepare_survey_data}* If performing an OOS validation of 
#' model performance, you should filter this dataset for the years "held back" 
#' from your model fit.
#' @param areas `sf` shapefile for specific country/region. Only required 
#' if `survey_circumcision_test` has records for area levels higher (i.e. more
#' granular) than `area_lev`, in which case they must be reassigned to their 
#' `parent_area_id` at `area_lev`, Default = NULL. 
#' @param area_lev Area level you wish to aggregate to when performing posterior
#' predictive comparisons with survey estimates.
#' @param type Decides type of circumcision coverage to perform PPC on, must 
#' be one of "MC", "MMC", or "TMC", Default = "MMC"
#' @param age_groups Age groups to aggregate by, Default:
#' c("0-4",   "5-9",   "10-14", "15-19", "20-24", "25-29",
#' "30-34", "35-39", "40-44", "45-49", "50-54", "54-59")
#' @param CI_range CI interval about which you want to compare empirical and
#' posterior predictive estimates for left out surveys, 
#' Default = c(0.5, 0.8, 0.95)
#' @param N Number of samples to generate, Default: 1000
#' @param seed Random seed used for taking binomial sample from posterior 
#' predictive distribution.
#' @return \code{data.frame} with samples aggregated by \code{aggr_cols} and
#' weighted by population.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname threemc_ppc
#' @export
threemc_ppc <- function(fit,
                        out,
                        survey_circumcision_test, 
                        areas = NULL,
                        area_lev = 1,
                        type = "MMC",
                        age_groups = c(
                          # five-year age groups
                          "0-4", "5-9", "10-14", "15-19",
                          "20-24", "25-29", "30-34", "35-39",
                          "40-44", "45-49", "50-54", "54-59" 
                        ),
                        CI_range = c(0.5, 0.8, 0.95),
                        N = 1000,
                        seed = 123) {
  
  stopifnot(type %chin% c("MC", "MMC", "TMC"))
  stopifnot("sample" %in% names(fit))
  
  # global bindings for data.table non-standard evaluation
  samp_colnames <- age <- age_group <- area_id <- area_level <- circ_status <- 
    indweight <- predicted <- sim_prop <- simulated <- year <- . <- 
    ..samp_colnames <- NULL
  
  # initialise data.table objects 
  out <- data.table::copy(data.table::setDT(out))
  survey_circumcision_test <- data.table::copy(
    data.table::setDT(survey_circumcision_test)
  )

  # filter results to specified or modelled area level
  if (!"area_level" %chin% names(out)) {
    out$area_level <- ifelse(
      nchar(out$area_id) == 3, 
      0, 
      substr(out$area_id, 4, 4)
    ) 
  }
  max_area_lev <- max(out$area_level)
  if (area_lev > max_area_lev) {
    message(
      "Specified area level more granular than data provided, area_lev set to ",
      max_area_lev
    )
    area_lev <- max_area_lev
  }
  
  # remove unneeded columns from survey_circumcision_test
  unneeded_cols <- c(
      "individual_id", "cluster_id", "household",
      "line", "sex", "dob_cmc", "interview_cmc"
  )
  unneeded_cols <- unneeded_cols[
    unneeded_cols %chin% names(survey_circumcision_test)
  ]
  survey_circumcision_test[, c(unneeded_cols) := NULL]
  
  # add year (i.e. survey_year) if not present
  if ("survey_year" %chin% names(survey_circumcision_test)) {
    survey_circumcision_test$year <- survey_circumcision_test$survey_year
  } else if (!"year" %chin% names(survey_circumcision_test)) {
    survey_circumcision_test$year <- as.numeric(
      substr(survey_circumcision_test$survey_id, 4, 7)
    )
  }


  #### Preprocess survey data and samples ####

  # Take sample matrix rows that are kept in .data from original skeleton data
  if ("n" %in% names(out)) {
    fit$sample <- lapply(fit$sample, function(x) x[out$n, ])
    out$n <- NULL
  } 
  # throw error if number of rows in results does not equal sample number
  stopifnot(nrow(out) == nrow(fit$sample$haz))
  
  # pull N "cum_inc" (only doing coverage) samples from fit
  samples <- fit$sample[
    grepl("cum_inc", names(fit$sample))
  ]
  # ensure names for MC columns in fit have the suffix "_mc"
  samples <- append_mc_name(samples)
  # remove "cum_inc_" and capitalise
  names(samples) <- paste(
    toupper(stringr::str_replace_all(names(samples), "cum_inc_", "")),
    "coverage"
  )
  rm(fit)

  # take samples of the correct "type"
  samples <- samples[[which(names(samples) == paste(type, "coverage"))]]

  # only take sample rows corresponding to subsetted out, if appropriate
  if ("n" %chin% names(out)) {
    samples <- samples[c(out$n), ]
    out <- out[, c("n") := NULL]
  }
  
  # if type == "MC", take both MMC and TMC as "non-missing" circumcision
  check_types <- switch(
    type, 
    "MC"  = c("MMC", "TMC"),
    type
  )
  # first remove records with missing circ_status
  survey_circumcision_test <- survey_circumcision_test[!is.na(circ_status)]
  # next, parse circumcision type from circ_who and circ_where
  survey_circumcision_test <- find_circ_type(survey_circumcision_test)
  # change type not matching the specified argument type to Missing
  survey_circumcision_test[,
    type := ifelse(
      type %chin% check_types & circ_status == 1,
      paste(type, "coverage"),
      "Missing"
    )
  ]
  # change missing circumcisions to non-circumcisions, for weighted mean
  # don't do so for MC, as we can have missing type circumcisions here!
  if (!type == "MC") {
    survey_circumcision_test[,
      circ_status := ifelse(
        type == "Missing",
        0,
        circ_status
      )
    ]  
  # convert all type to MC
  } else {
    survey_circumcision_test[,
      type := ifelse(
        circ_status == 1,
        "MC coverage",
        type
      )
    ]  
  }
  
  
  #### Join Samples with Out and Aggregate to area_lev ####
  
  # take random sample of cols in samples if required
  if (N != ncol(samples)) {
    set.seed(seed)
    choose_cols <- sample(seq_len(ncol(samples)), size = N)
    samples <- samples[, choose_cols]
  }

  # join coverage col of interest with samples
  out_types <- out[, c("area_id", "area_level", "year", "age", "population")]
  
  n <- length(out_types)
  # much faster than dplyr::bind_cols
  out_types <- cbind(out_types, samples)
  out_types$type <- paste(type, "coverage")
  rm(out)
  
  samp_colnames <- sprintf("samp_%0*d", nchar(N), seq_len(N))
  
  # change col names to work in aggregate_sample_age_group
  names(out_types)[grepl("V", names(out_types))] <- samp_colnames
  
  out_types <- dplyr::relocate(
    out_types, type, .before = dplyr::contains("samp")
  )
  
  # reassign to desired `area_lev` (uses `<=` as we can't disaggregate!)
  if (!all(out_types$area_level <= area_lev)) {
    # must provide areas for this
    if (is.null(areas)) stop("areas required for reassigning area level")

    # reassign area levels in out to `area_lev`
    out_types <- reassign_survey_level(out_types, areas, area_lev)

    # aggregate num_cols (i.e. incidence, rate, coverage) in out by aggr_cols
    out_cols <- names(out_types)
    out_cols <- out_cols[!out_cols %chin% c("space", "population")]
    aggr_cols <- out_cols[seq_len(which(out_cols == "type"))]
    # out_types <- dplyr::as_tibble(aggregate_sample(
    out_types <- aggregate_sample(
      .data     = out_types,
      aggr_cols = aggr_cols,
      num_cols = out_cols[!out_cols %chin% c(aggr_cols)]
    )
  }
  # filter for test year(s) and modelled area level
  out_types <- out_types[
    year %in% survey_circumcision_test$year & area_level == area_lev
  ]

  
  #### Prepare Survey Points & Join with Samples ####

  # reassign test data area level to area_lev, if any are more granular
  if (any(survey_circumcision_test$area_level > area_lev)) {
    if (is.null(areas)) stop("areas required for reassigning area level")
    
    survey_circumcision_test <- reassign_survey_level(
        survey_circumcision = survey_circumcision_test,
        areas               = areas,
        area_lev            = area_lev
    )
  }

  # filter to modelled area_lev and ages
  survey_estimate_prep <- survey_circumcision_test[
    area_level == area_lev & age < max(out_types$age)
  ]

  # join with samples
  out_types <- out_types[
    , 
    c("area_id", "year", "age", samp_colnames), 
    with = FALSE
  ]
  survey_estimate_ppd <- data.table::merge.data.table(
    survey_estimate_prep, out_types, all.x = TRUE
  )
  
  # stop (or give message) for unusable survey_estimate_ppd/NAs in sample cols
  if (nrow(survey_estimate_ppd) == 0) {
    stop(
      "No matches between test data (survey_circumcision_test) and model",
      " results, (out), check inputs"
    )
  }
  samp_cols <- grepl("samp_", names(survey_estimate_ppd))
  if (any(is.na(survey_estimate_ppd[, samp_cols, with = FALSE]))) {
    message(
      "Some NAs when matching PPD samples and test survey data, check ",
      "that ages and years match in both"
    )
  }
  rm(survey_circumcision_test, survey_estimate_prep, out_types)
  
  
  #### Binomial Sample from PPD ####

  # switch samp columns to long "predicted" column 
  survey_estimate_ppd_long <- data.table::melt(
    survey_estimate_ppd, 
    measure       = patterns("^samp"), 
    variable.name = "sample",
    value.name    = "predicted"
  )
  
  set.seed(seed)
  survey_estimate_ppd_long[,
    predicted := ifelse(predicted > 1, 1, predicted)
  ][, 
    simulated := stats::rbinom(
     nrow(survey_estimate_ppd_long), 1L, prob = predicted
    )
  ]
  gc()
  
  # give message about NAs
  if (any(is.na(survey_estimate_ppd_long$simulated))) {

    na_sims <- survey_estimate_ppd_long[
      is.na(simulated)
    ][, 
      .N, 
      by = c("survey_id", "area_id", "year", "type")
    ][
      order(year, area_id, type)
    ]

    message("The following have had NAs removed:")
    message(
      paste0(utils::capture.output(data.frame(na_sims)), collapse = "\n")
    )
    rm(na_sims)
   
    survey_estimate_ppd_long <- survey_estimate_ppd_long[!is.na(simulated)]
  }
  
  
  #### Group by Age Group ####
  
  # create df to match five year age groups to single ages
  age_group_df <- data.table::rbindlist(lapply(
    age_groups,
    match_age_group_to_ages,
    max_age = max(survey_estimate_ppd_long$age)
  ))

  # join in age groups matching each single age
  survey_estimate_age_group <- data.table::merge.data.table(
    survey_estimate_ppd_long, 
    age_group_df, 
    all.x = TRUE, 
    allow.cartesian = TRUE
  )[
    # remove missing age groups and unknown circ_status
    !is.na(circ_status) & !is.na(age_group)
  ]
  
  # calc survey mean separately (will be repeated for each survey otherwise)
  survey_estimate_age_group_mean <- survey_estimate_age_group[
    ,
    .(mean = stats::weighted.mean(x = circ_status, w = indweight)),
    by = c("area_id", "year", "age_group")
  ]
  
  # calculate simulated proportions
  survey_estimate_age_group <- survey_estimate_age_group[
    ,
    .(sim_prop = stats::weighted.mean(simulated, indweight)),
    by = c("area_id", "year", "age_group", "sample")
  ][
    ,
    type := paste(type, "coverage")
  ]
  # add sample mean to simulated proportions
  survey_estimate_age_group <- data.table::merge.data.table(
    survey_estimate_age_group,
    survey_estimate_age_group_mean, 
    all.x = TRUE
  )
  survey_estimate_age_group <- dplyr::relocate(
    survey_estimate_age_group, type, .before = "sim_prop"
  )
  
  rm(age_group_df, survey_estimate_age_group_mean)
  
  # give warning about missing age groups
  missing_age_groups <- age_groups[
    !age_groups %chin% survey_estimate_age_group$age_group
  ]
  if (length(missing_age_groups) > 0) {
    message(paste0(
      "The following `age_groups` are missing in `survey_estimate`:\n",
      paste(missing_age_groups, collapse = ", ")
    ))
  }
  
  # remove any NAs (because age groups and/or years missing for some areas)
  if (any(is.na(survey_estimate_age_group$mean))) {
    na_surveys <- survey_estimate_age_group[
      is.na(mean)
    ][, 
      .N, 
      by = c("area_id", "year", "age_group", "type")
    ][
      order(year, age_group, area_id)
    ]
    
    message("The following have had NAs removed:")
    message(
      paste0(utils::capture.output(data.frame(na_surveys)), collapse = "\n")
    )
    rm(na_surveys)
    
    survey_estimate_age_group <- survey_estimate_age_group[!is.na(mean)]
  }
  
  # find quantiles for age group PPDs
  ppd_quantiles <- survey_estimate_age_group[
    ,
    .(ppd_mean = mean(sim_prop),
      ppd_0.025  = stats::quantile(sim_prop, 0.025),
      ppd_0.100  = stats::quantile(sim_prop, 0.100),
      ppd_0.250  = stats::quantile(sim_prop, 0.250),
      ppd_median = stats::quantile(sim_prop, 0.500),
      ppd_0.750  = stats::quantile(sim_prop, 0.750),
      ppd_0.900  = stats::quantile(sim_prop, 0.900),
      ppd_0.975  = stats::quantile(sim_prop, 0.975)),
    by = c("area_id", "year", "age_group", "type")
  ]
  
  rm(survey_estimate_ppd_long)
  
  
  #### Posterior Predictive Check for Prevalence Estimations ####

  # fun to calculate where in model sample distribution empirical values are
  quant_pos_sum <- function(y, x) if (y < x) 0 else 1

  # convert back to wide format
  survey_estimate_ppd_wide <- data.table::dcast(
    survey_estimate_age_group, ... ~ sample, value.var = "sim_prop"
  )

  # find position of mean estimate from surveys amongst predictions
  survey_estimate_ppd_dist <- survey_estimate_ppd_wide[
    , 
    lapply(.SD, function(x) quant_pos_sum(mean, x)),
    by = c("area_id", "year", "age_group", "type"),
    .SDcols = samp_colnames
  ]
  survey_estimate_ppd_dist$quant_pos <- apply(
    survey_estimate_ppd_dist[, c(samp_colnames), with = FALSE], 1, sum
  )
  # remove sample columns
  survey_estimate_ppd_dist <- survey_estimate_ppd_dist[
    , c(samp_colnames) := NULL
  ]
  
  # add quant pos columns to dataframe with survey obs and PPD samples
  survey_estimate_ppd <- data.table::merge.data.table(
    survey_estimate_ppd_dist, survey_estimate_ppd_wide, all.x = TRUE
  )
  
  # calculate position of oos obs within ordered PPD (i.e. estimate of hist)
  oos_within_ppd <- function(x, CI_range) {
    sum(
      x >= 0.5 * (1 - CI_range) * N &
        x <= N - 0.5 * (1 - CI_range) * N
    ) / length(x)
  }
  CI_range <- sort(CI_range)
  oos_within_ppd_percent <- vapply(CI_range, function(x) {
    oos_within_ppd(survey_estimate_ppd$quant_pos, x)
  }, numeric(1))
  names(oos_within_ppd_percent) <- CI_range

  # assign to message so we don't print to console with [[i]] separators
  message <- lapply(seq_along(CI_range), function(i) {
    print(paste0(
      "Percentage of survey points which fall within posterior predictive",
      " distribution at ",
      CI_range[[i]] * 100,
      "%",
      " CI: ",
      round(oos_within_ppd_percent[[i]][1] * 100, 2),
      "%"
    ))
  })

  summary_stats <- list(
    "oos_observations_within_PPD_CI" = oos_within_ppd_percent
  )
  
  # join PPD quantiles with quant_pos and samples
  survey_estimate_ppd <- data.table::merge.data.table(
    survey_estimate_ppd, ppd_quantiles, all.x = TRUE
  )
  
  rm(
    survey_estimate_age_group,
    survey_estimate_ppd_wide,
    survey_estimate_ppd_dist,
    ppd_quantiles
  )
  gc()
  

  #### Calculate ELPD, CRPS & Error Stats ####

  # compute summary stats for comparison with other models, if specified
  # Actual OOS survey obs vs posterior predictive distribution for each
  actual <- survey_estimate_ppd$mean
  predictions <- survey_estimate_ppd[, ..samp_colnames]
  
  # calculate elpd, crps, mean error stats
  actual_pred_comparison <- compare_predictions(actual, predictions)
  
  # Add pointwise error stats to return df
  survey_estimate_ppd <- cbind(
    survey_estimate_ppd, actual_pred_comparison$summary_stats_df
  )
  
  # add mean error stats to return list
  summary_stats <- c(
    summary_stats,
    actual_pred_comparison$summary_stats
  )
  
  # move sample columns to the end (also moves summary cols to the start)
  survey_estimate_ppd <- survey_estimate_ppd %>%
    dplyr::relocate(dplyr::contains("samp"), .after = dplyr::everything())

  # Return results
  return(list(
    "ppc_df"        = survey_estimate_ppd,
    "summary_stats" = summary_stats
  ))
}

#' @title Define circumcision type
#' @description Using `circ_who` and `circ_where`, determines survey 
#' type. 
#' @param actual `vector` of observed values. 
#' @param predictions `dataframe/tibble/data.table` or `matrix` with 
#' a row for each value in `actual`, and a column for each prediction. 
#' @return List of: 
#' \itemize{
#'  \item{"summary_stats_df"}{pointwise prediction stats as a `data.frame`}
#'  \item{"summary_stats"}{mean/summed overall prediction stats}
#' }
#' @export
#' @rdname reassign_surey_level
#' @keywords internal
compare_predictions <- function(actual, predictions) {
  
  # convert to matrix (required for error calculations)
  predictions <- as.matrix(predictions)
  
  # calculate ELPD
  elpd <- loo::elpd(t(predictions))
  
  # calculate CRPS
  crps <- scoringutils::crps_sample(
    true_values = actual,
    predictions = as.matrix(predictions)
  )
  
  # calculate error stats (MAE, MSE & RMSE)
  
  # errors for each sample from PPD for each observation
  errors <- actual - predictions # errors
  ae <- abs(errors) # absolute errors
  se <- errors ^ 2 # squared error
  
  # mean errors for each observations
  mae_vec <- apply(ae, 1, mean)
  mse_vec <- apply(se, 1, mean)
  rmse_vec <- sqrt(mse_vec)
  
  # mean errors overall
  mae <- mean(mae_vec)
  mse <- mean(mse_vec)
  rmse <- sqrt(mse)
  
  return(list(
    # pointwise summaries
    "summary_stats_df" = data.frame(
      "elpd" = elpd$pointwise[, 1], # pointwise elpd
      "crps" = crps,
      "mae"  = mae_vec,
      "mse"  = mse_vec,
      "rmse" = rmse_vec
    ), 
    # mean summaries
    "summary_stats" = list(
      "elpd" = elpd,
      "crps" = sum(crps),
      "mae"  = mae,
      "mse"  = mse,
      "rmse" = rmse
    )
  ))
}
