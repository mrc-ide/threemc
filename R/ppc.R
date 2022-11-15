#' @title Posterior Predictive Distribution and checks on OOS survey
#' @description Aggregate specified `numeric` columns by population-weighted
#' age groups (rather than single year ages), split by specified categories.
#' @param fit Fit object returned by \link[naomi]{sample_tmb}, which includes,
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
#' @param age_groups Age groups to aggregate by, Default:
#' c("0-4",   "5-9",   "10-14", "15-19", "20-24", "25-29",
#' "30-34", "35-39", "40-44", "45-49", "50-54", "54-59")
#' @param CI_range CI interval about which you want to compare empirical and
#' posterior predictive estimates for left out surveys, 
#' Default = c(0.5, 0.8, 0.95)
#' @param N Number of samples to generate, Default: 1000
#' @param compare_stats Set to TRUE if you wish to compute comparative
#' statistics (specifically, ELPD and CRPS) to compare with alternative models,
#' Default: TRUE
#' @return \code{data.frame} with samples aggregated by \code{aggr_cols} and
#' weighted by population.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @rdname threemc_oos_ppc
#' @export
threemc_oos_ppc <- function(fit,
                            out,
                            survey_circumcision_test, 
                            areas = NULL,
                            area_lev = 1,
                            age_groups = c(
                              # five-year age groups
                              "0-4", "5-9", "10-14", "15-19",
                              "20-24", "25-29", "30-34", "35-39",
                              "40-44", "45-49", "50-54", "54-59" 
                            ),
                            CI_range = c(0.5, 0.8, 0.95),
                            N = 1000,
                            compare_stats = TRUE) {


  # check that all age groups are in survey estimates
  # TODO: Must check this later with test data
  # missing_age_groups <- age_groups[!age_groups %chin% survey_estimate$age_group]
  # if (length(missing_age_groups) > 0) {
  #   message(paste0(
  #     "The following `age_groups` are missing in `survey_estimate`:\n",
  #     paste(missing_age_groups, collapse = ", ")
  #   ))
  # }


  # remove unneeded columns from survey_circumcision_test
  unneeded_cols <- c(
      "individual_id", "cluster_id", "household",
      "line", "sex", "dob_cmc", "interview_cmc"
  )
  survey_circumcision_test <- survey_circumcision_test %>%
      dplyr::select(-dplyr::any_of(unneeded_cols))
  
  # add year (i.e. survey_year) if not present
  if ("survey_year" %in% names(survey_circumcision_test)) {
    survey_circumcision_test$year <- survey_circumcision_test$survey_year
  } else if (!"year" %in% names(survey_circumcision_test)) {
    survey_circumcision_test <- survey_circumcision_test %>% 
      dplyr::mutate(year = as.numeric(substr(.data$survey_id, 4, 7)))
  }


  #### Join Samples with Results ####

  # filter results to specified or modelled area level
  max_area_lev <- max(out$area_level)
  if (area_lev > max_area_lev) {
    message(
      "Specified area lev more granular than data provided, set to ",
      max_area_lev
    )
    area_lev <- max_area_lev
  }

  # filter to area_lev
  out <- dplyr::filter(out, .data$area_level == area_lev)
  # join in populations
  # dplyr::left_join(
  #   dplyr::select(
  #     populations,
  #     -c(dplyr::matches("area_name"), dplyr::matches("area_level"))
  #   )
  # )

  # pull sample name corresponding with type column values in results
  # sample_colnames <- switch(type,
  #   "coverage"    = "cum_inc",
  #   "incidence"   = "inc",
  #   "probability" = "haz"
  # )
  # if (is.null(sample_colnames)) {
  #   stop(
  #     "Please choose a valid type
  #       (one of 'probability', 'incidence', 'prevalence'"
  #   )
  # }

  # pull N "cum_inc" samples for each desired output from fit
  samples <- fit$sample[
    # grepl(paste0("^", sample_colnames), names(fit$sample))
    grepl("cum_inc", names(fit$sample))
  ]
  # ensure names for MC columns in fit have the suffix "_mc"
  samples <- append_mc_name(samples)
  # remove "cum_inc_" and capitalise
  names(samples) <- paste(
    toupper(stringr::str_replace_all(names(samples), "cum_inc_", "")),
    "coverage"
  )

  # TEMP: take only medical
  samples <- samples[names(samples) == "MMC coverage"]
  survey_circumcision_test <- find_circ_type(survey_circumcision_test) %>%
      # filter(type %in% c("MMC", "Missing")) %>%
      # mutate(type = "MMC Coverage")
      # Any non-medical circumcision is treated as a "competing risk"
      dplyr::mutate(type = ifelse(.data$type == "MMC", "MMC Coverage", "Missing"))

  #### Split out between types ####

  # if (is.null(removed_years)) {
  #   message(
  #     "No `removed_years` specified, performing non-OOS PPC using all surveys"
  #   )
  #   removed_years <- unique(out$year)
  # }

  # join coverage col of interest with samples
  # NOTE: Only doing for medical now, so can change
  out_types <- lapply(names(samples), function(x) {
    out_spec <- dplyr::select(out, .data$area_id:.data$population)
    n <- length(out_spec)
    out_spec[, (n + 1):(n + N)] <- samples[[x]]
    out_spec <- dplyr::mutate(out_spec, indicator = x)
  }) %>%
    dplyr::bind_rows() %>%
    # filter for test year(s) and modelled area level
    # dplyr::filter(.data$year %in% removed_years, .data$area_level == area_lev)
    dplyr::filter(
      .data$year %in% survey_circumcision_test$year,
      .data$area_level == area_lev
    )

  ## out_types <- dplyr::select(out, .data$area_id:.data$population)
  ## n <- length(out_types)
  ## out_types[, (n + 1):(n + N)] <- samples
  ## out_types <- dplyr::mutate(out_types, indicator = "MMC coverage") %>%
  ##     # filter for test year(s) and modelled area level
  ##     dplyr::filter(
  ##         .data$year %in% survey_circumcision_test$year,
  ##         .data$area_level == area_lev
  ##     )

  # change col names to work in aggregate_sample_age_group
  names(out_types)[grepl("V", names(out_types))] <- paste0("samp_", 1:N)
  # out_types <- dplyr::rename(out_types, type = .data$indicator)
  out_types <- out_types %>%
    dplyr::relocate(type = .data$indicator, .before = dplyr::contains("samp"))


  #### Aggregate to Age Groups ####

  # out_types_agegroup <- aggregate_sample_age_group(
  #   out_types,
  #   aggr_cols = c(
  #     "area_id", "area_name", "area_level", "year", "type"
  #   ),
  #   num_cols = paste0("samp_", seq_len(N)),
  #   age_groups = age_groups
  # ) %>%
  #   # rename to match survey points df
  #   dplyr::rename(indicator = .data$type) %>%
  #   dplyr::relocate(
  #     dplyr::contains("samp_"),
  #     .after = dplyr::everything()
  #   ) %>%
  #   # filter for at least area_level 1
  #   dplyr::filter(.data$area_level == min(area_lev, 1))
  #
  # # check for NAs in posterior predictive distributions
  # stopifnot(!all(is.na(
  #   as.matrix(dplyr::select(out_types_agegroup, dplyr::contains("samp")))
  # )))


  #### Prepare Survey Points & Join with Samples ####

  # survey_estimate_prep <- survey_estimate %>%
  #   # filter for OOS year(s) and modelled area_level
  #   dplyr::filter(
  #     .data$year %in% removed_years,
  #     .data$area_level == min(area_lev, 1),
  #     .data$age_group %chin% out_types_agegroup$age_group,
  #     .data$mean != 0
  #   ) %>%
  #   # ignore survey_id, merging for the same year
  #   dplyr::group_by(
  #     .data$area_id,
  #     .data$year,
  #     .data$age_group,
  #     .data$indicator
  #   ) %>%
  #   dplyr::summarise(dplyr::across(
  #     # dplyr::all_of(c("mean", "upper", "lower")),
  #     dplyr::all_of(c("mean", "upper", "lower", "n_eff_kish")),
  #     sum,
  #     na.rm = TRUE
  #   ), .groups = "drop") %>%
  #   # floor kish effective sample size
  #   dplyr::mutate(n_eff_kish = floor(.data$n_eff_kish))

  # reassign test data area level to area_lev, if any are more granular
  if (any(survey_circumcision_test$area_level > area_lev)) {
    # must provide areas for this
    if (is.null(areas)) stop("areas required for reassigning area level")
    # remove unwanted geometry information, if present
    if (inherits(areas, "sf")) areas <- sf::st_drop_geometry(areas)
    survey_circumcision_test <- reassign_survey_level(
        survey_circumcision = survey_circumcision_test,
        areas               = dplyr::select(areas, dplyr::contains("area")),
        area_lev            = area_lev
    )
  }

  # filter to modelled area_lev and ages
  survey_estimate_prep <- survey_circumcision_test %>%
    dplyr::filter(
      .data$area_level == area_lev, 
      .data$age < max(out_types$age)
    )

  # join with samples
  survey_estimate_ppd <- survey_estimate_prep %>%
    dplyr::left_join(
      # out_types_agegroup
      out_types %>%
        dplyr::select(.data$area_id, .data$year, .data$age, dplyr::starts_with("samp"))
    )

  # stop (or give message) for unusable survey_estimate_ppd/NAs in sample cols
  if (nrow(survey_estimate_ppd) == 0) {
    stop(
      "No matches between test data (survey_circumcision_test) and model results,",
      " (out), check inputs"
    )
  }
  samp_cols <- grepl("samp_", names(survey_estimate_ppd))
  if (any(is.na(survey_estimate_ppd[, samp_cols]))) {
    message(
      "Some NAs when matching PPD samples and test survey data, check ",
      "that ages and years match in both"
    )
  }


  #### Binomial Sample from PPD ####

  # pull out columns with samples and n_eff_kish
  # samp_df <- survey_estimate_ppd %>%
  #   dplyr::select(.data$n_eff_kish, dplyr::contains("samp_"))
  #
  # # take N samples from n_eff_kish trials with success probability "samp_x"
  # samples_binom <- data.frame(t(apply(
  #   samp_df, 1, function(x) {
  #     stats::rbinom(N, x["n_eff_kish"], x[names(x) != "n_eff_kish"])
  #   }
  # )))
  # names(samples_binom) <- paste0("samp_", seq_len(N))
  #
  # survey_estimate_ppd_dist <- survey_estimate_ppd %>%
  #   # replace coverage sample columns with binomial samples
  #   dplyr::select(-dplyr::contains("samp_")) %>%
  #   dplyr::bind_cols(samples_binom) %>%
  #   # divide binomial samples by n_eff_kish; compare these to survey estimates
  #   dplyr::mutate(
  #     dplyr::across(dplyr::contains("samp_"), ~ .x / .data$n_eff_kish)
  #   ) %>%
  #   dplyr::ungroup() %>%

  # switch samp columns to long "predicted" column 
  # Could probably do something less memory intensive here!
  set.seed(123)
  survey_estimate_ppd_long <- survey_estimate_ppd %>%
    tidyr::pivot_longer(
      dplyr::starts_with("samp"),
      names_to  = "sample",
      values_to = "predicted"
    ) %>%
    # sample at individual level from binomial dist with success prob from PPD
    dplyr::mutate(
      simulated = stats::rbinom(dplyr::n(), 1, prob = .data$predicted)
    )

  
  #### Group by Age Group ####

  # create df to match five year age groups to single ages
  age_group_df <- dplyr::bind_rows(lapply(
    age_groups,
    match_age_group_to_ages,
    max_age = max(survey_estimate_ppd_long$age)
  ))

  survey_estimate_age_group <- survey_estimate_ppd_long %>%
    # join in age groups matching each single age
    dplyr::left_join(age_group_df, by = "age") %>%
    # remove missing age groups and unknown circ_status
    dplyr::filter(!is.na(.data$circ_status), !is.na(.data$age_group)) %>%
    # aggregate single ages to age groups
    dplyr::group_by(.data$area_id, .data$year, .data$age_group, .data$sample) %>%
    # calculate weighted means for observed and simulated proportions
    dplyr::summarise(
      # obs_prop = stats::weighted.mean(circ_status, indweight), # replace
      mean     = stats::weighted.mean(.data$circ_status, .data$indweight), # better to do this with wide format, recalculating the same value each time here!
      sim_prop = stats::weighted.mean(.data$simulated, .data$indweight),
      .groups = "drop"
    )

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
  
  # find quantiles for age group PPDs
  # TODO: Can do a better across here!
  ppd_quantiles <- survey_estimate_age_group %>%
    # survey_estimate_age_group %>%
    dplyr::group_by(.data$area_id, .data$year, .data$age_group) %>%
    dplyr::summarise(
      # mean       = mean(mean),
      ppd_mean   = mean(.data$sim_prop),
      ppd_0.025  = stats::quantile(.data$sim_prop, 0.025),
      ppd_0.10   = stats::quantile(.data$sim_prop, 0.1),
      ppd_0.25   = stats::quantile(.data$sim_prop, 0.25),
      ppd_median = stats::quantile(.data$sim_prop, 0.5),
      ppd_0.75   = stats::quantile(.data$sim_prop, 0.75),
      ppd_0.90   = stats::quantile(.data$sim_prop, 0.9),
      ppd_0.975  = stats::quantile(.data$sim_prop, 0.975)
    )
  

  #### Posterior Predictive Check for Prevalence Estimations ####

  # find middle samp column
  # samp_cols <- as.numeric(gsub(
  #   ".*?([0-9]+).*", "\\1",
  #   names(survey_estimate_ppd)[grepl("samp", names(survey_estimate_ppd))]
  # ))
  # mid_samp <- samp_cols[length(samp_cols) / 2]

  # fun to calculate where in model sample distribution empirical values are
  quant_pos_sum <- function(y, x) if (y < x) 0 else 1

  survey_estimate_ppd_wide <- survey_estimate_age_group %>%
    # convert back to wide format
    tidyr::pivot_wider(names_from = "sample", values_from = "sim_prop")

  survey_estimate_ppd_dist <- survey_estimate_ppd_wide %>%
    dplyr::group_by(dplyr::across(.data$area_id:.data$age_group)) %>%
    dplyr::summarise(
      # find position of mean estimate from surveys amongst predictions
      quant_pos = sum(
        dplyr::across(dplyr::starts_with("samp_"), ~ quant_pos_sum(
          mean, .x
        ))
      ),
      .groups = "rowwise"
    ) %>%
    dplyr::ungroup()

  # add quant pos columns to dataframe with survey obs and PPD samples
  survey_estimate_ppd <- survey_estimate_ppd_dist %>%
    dplyr::left_join(survey_estimate_ppd_wide, by = c("area_id", "year", "age_group"))

  # calculate position of oos obs within ordered PPD (i.e. estimate of hist)
  oos_within_ppd <- function(x, CI_range) {
    sum(
      x >= 0.5 * (1 - CI_range) * N &
        x <= N - 0.5 * (1 - CI_range) * N
    ) / length(x)
  }
  CI_range <- sort(CI_range)
  oos_within_ppd_percent <- vapply(CI_range, function(x) {
    # oos_within_ppd(survey_estimate_ppd$quant_pos_final, x)
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
  survey_estimate_ppd <- survey_estimate_ppd %>%
    dplyr::left_join(ppd_quantiles) 
  

  #### Calculate ELPD, CRPS & Error Stats ####

  # compute summary stats for comparison with other models, if specified
  if (compare_stats == TRUE) {

    # Actual OOS survey obs vs posterior predictive distribution for each
    # actual <- survey_estimate_ppd$mean
    actual <- survey_estimate_ppd$mean
    predictions <- dplyr::select(survey_estimate_ppd, dplyr::contains("samp_"))

    # calculate ELPD
    elpd <- loo::elpd(t(predictions))
    
    # calculate CRPS
    crps <- scoringutils::crps_sample(
      true_values = actual,
      predictions = as.matrix(predictions)
    )
    
    # add pointwise elpd and crps to output df
    survey_estimate_ppd$elpd <- elpd$pointwise[, 1]
    survey_estimate_ppd$crps <- crps

    # calculate error stats (MAE, MSE & RMSE)

    # errors for each sample from PPD for each observation
    errors <- actual - as.matrix(predictions) # errors
    ae <- abs(errors) # absolute errors
    se <- errors^2 # squared error

    # mean errors for each observations
    mae_vec <- apply(ae, 1, mean)
    mse_vec <- apply(se, 1, mean)
    rmse_vec <- sqrt(mse_vec)

    # mean errors overall
    mae <- mean(mae_vec)
    mse <- mean(mse_vec)
    rmse <- sqrt(mse)

    # Add error stats to return df
    survey_estimate_ppd <- survey_estimate_ppd %>%
      dplyr::mutate(
        mae = mae_vec,
        mse = mse_vec,
        rmse = rmse_vec
      )

    summary_stats <- c(
      summary_stats,
      list(
        "elpd" = elpd,
        "crps" = crps,
        "mae"  = mae,
        "mse"  = mse,
        "rmse" = rmse
      )
    )
  }

  survey_estimate_ppd <- survey_estimate_ppd %>%
    dplyr::relocate(dplyr::contains("samp"), .after = dplyr::everything())

  # Return results
  return(list(
    "ppc_df"        = survey_estimate_ppd,
    "summary_stats" = summary_stats
  ))
}
