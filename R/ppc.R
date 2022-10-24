#' @title Posterior Predictive Distribution and checks on OOS survey
#' @description Aggregate specified `numeric` columns by population-weighted
#' age groups (rather than single year ages), split by specified categories.
#' @param fit Fit object returned by \link[naomi]{sample_tmb}, which includes,
#' among other things, the optimised parameters and subsequent sample for our
#' TMB model.
#' @param out Unaggregated results of model fitting outputted by
#' \link[threemc]{compute_quantiles}.
#' @param populations Single age male population counts by space and time.
#' @param survey_estimate Circumcision estimates based on empirical observations
#' from surveys.
#' @param removed_years Years removed from dataset used to fit model
#' (i.e. "withheld" surveys from OOS model fitting).
#' @param type The desired circumcision estimate. Can be one of
#' "probability", "incidence" or "coverage".
#' @param area_lev Area level you wish to aggregate to when performing posterior
#' predictive comparisons with survey estimates.
#' @param age_groups Age groups to aggregate by, Default:
#' c("0-4",   "5-9",   "10-14", "15-19", "20-24", "25-29",
#' "30-34", "35-39", "40-44", "45-49", "50-54", "54-59",
#' "0+",    "10+",   "15+",   "15-24", "10-24", 15-29",
#' "10-29", "15-39", "10-39", "15-49", "10-49")
#' @param CI_range CI interval about which you want to compare empirical and
#' posterior predictive estimates for left out surveys.
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
                            populations,
                            survey_estimate,
                            removed_years,
                            type = "coverage",
                            area_lev = 1,
                            age_groups = c(
                              # five-year age groups
                              "0-4",   "5-9",   "10-14", "15-19", 
                              "20-24", "25-29", "30-34", "35-39", 
                              "40-44", "45-49", "50-54", "54-59",
                              # age groups with only minimum cut-off
                              "0+", "10+", "15+",
                              # other, wider age groups of interest
                              "10-24", "15-24", "10-29", "15-29",
                              "10-39", "15-39", "10-49", "15-49"
                            ),
                            CI_range = 0.95,
                            N = 1000,
                            compare_stats = TRUE) {

  
  # check that all age groups are in survey estimates
  missing_age_groups <- age_groups[!age_groups %chin% survey_estimate$age_group]
  if (length(missing_age_groups) > 0) {
    message(paste0(
      "The following `age_groups` are missing in `survey_estimate`:\n",
      paste(missing_age_groups, collapse = ", ")
    ))
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
  out <- dplyr::filter(out, .data$area_level == area_lev) %>%
    # join in populations
    dplyr::left_join(
      dplyr::select(
        populations,
        -c(dplyr::matches("area_name"), dplyr::matches("area_level"))
      )
    )
  
  # pull sample name corresponding with type column values in results
  sample_colnames <- switch(type,
    "coverage"    = "cum_inc",
    "incidence"   = "inc",
    "probability" = "haz"
  )
  if (is.null(sample_colnames)) {
    stop(
      "Please choose a valid type
        (one of 'probability', 'incidence', 'prevalence'"
    )
  }

  # pull N "cum_inc" samples for each desired output from fit
  samples <- fit$sample[
    grepl(paste0("^", sample_colnames), names(fit$sample))
  ]
  # ensure names for MC columns in fit have the suffix "_mc"
  samples <- append_mc_name(samples)
  # remove "cum_inc_" and capitalise
  names(samples) <- paste(
    toupper(stringr::str_replace_all(names(samples), "cum_inc_", "")), type
  )


  #### Split out between types ####

  # join coverage col of interest with samples
  out_types <- lapply(names(samples), function(x) {
    out_spec <- dplyr::select(out, .data$area_id:.data$population)
    n <- length(out_spec)
    out_spec[, (n + 1):(n + N)] <- samples[[x]]
    out_spec <- dplyr::mutate(out_spec, indicator = x)
  }) %>%
    dplyr::bind_rows() %>%
    # only take years where surveys were removed, and modelled area level
    dplyr::filter(.data$year %in% removed_years, .data$area_level == area_lev)


  #### Aggregate to Age Groups ####

  # change col names to work in aggregate_sample_age_group
  names(out_types)[grepl("V", names(out_types))] <- paste0("samp_", 1:N)
  out_types <- dplyr::rename(out_types, type = .data$indicator)

  out_types_agegroup <- aggregate_sample_age_group(
    out_types,
    aggr_cols = c(
      "area_id", "area_name", "area_level", "year", "type"
    ),
    num_cols = paste0("samp_", seq_len(N)),
    age_groups = age_groups
  ) %>%
    # rename to match survey points df
    dplyr::rename(indicator = .data$type) %>%
    dplyr::relocate(
      dplyr::contains("samp_"),
      .after = dplyr::everything()
    ) %>%
    # filter for at least area_level 1
    dplyr::filter(.data$area_level == min(area_lev, 1))
  
  # check for NAs in posterior predictive distributions
  stopifnot(!all(is.na(
    as.matrix(dplyr::select(out_types_agegroup, dplyr::contains("samp")))
  )))


  #### Prepare Survey Points & Join with Samples ####

  survey_estimate_prep <- survey_estimate %>%
    # filter for OOS year(s) and modelled area_level
    dplyr::filter(
      .data$year %in% removed_years,
      .data$area_level == min(area_lev, 1),
      .data$age_group %chin% out_types_agegroup$age_group,
      .data$mean != 0
    ) %>%
    # ignore survey_id, merging for the same year
    dplyr::group_by(
      .data$area_id,
      .data$year,
      .data$age_group,
      .data$indicator
    ) %>%
    dplyr::summarise(dplyr::across(
      dplyr::all_of(c("mean", "upper", "lower")),
      sum, 
      na.rm = TRUE
    ), .groups = "drop")

  # join with samples
  survey_estimate_ppd <- dplyr::left_join(
    survey_estimate_prep, out_types_agegroup
  )
  

  #### Calculate Posterior Predictive Check for Prevalence Estimations ####

  # find middle samp column 
  samp_cols <- as.numeric(gsub(
    ".*?([0-9]+).*", "\\1", 
    names(survey_estimate_ppd)[grepl("samp", names(survey_estimate_ppd))]
  ))
  mid_samp <- samp_cols[length(samp_cols) / 2] 
  
  
  # func calculating where in model sample distribution empirical values are
  quant_pos_sum <- function(y, x) if (y < x) 0 else 1
  
  survey_estimate_ppd_dist <- survey_estimate_ppd %>%
    dplyr::relocate(
      dplyr::all_of(c("mean", "upper", "lower")), .before = .data$samp_1
    ) %>%
    dplyr::group_by(dplyr::across(.data$area_id:.data$area_level)) %>%
    dplyr::summarise(
      # find position of mean estimate from surveys amongst PPD
      quant_pos = sum(
        dplyr::across(dplyr::starts_with("samp_"), ~ quant_pos_sum(
          mean, .x
        ))
      ),
      # find corresponding uncertainty bounds (based on survey uncertainty)
      quant_pos_lower = sum(
        dplyr::across(dplyr::starts_with("samp_"), ~ quant_pos_sum(
          lower, .x
        ))
      ),
      quant_pos_upper = sum(
        dplyr::across(dplyr::starts_with("samp_"), ~ quant_pos_sum(
          upper, .x
        ))
      ),
      .groups = "rowwise"
    ) %>%
    # find optimum quant position (i.e. closest to middle of sample)
    tidyr::pivot_longer(dplyr::contains("quant_pos")) %>% 
    dplyr::mutate(diff_mid_samp = abs(.data$value - mid_samp)) %>% 
    dplyr::group_by(dplyr::across(.data$area_id:.data$area_level)) %>% 
    dplyr::mutate(
      quant_pos_final = .data$value[which.min(.data$diff_mid_samp)]
    ) %>%
    dplyr::select(-.data$diff_mid_samp) %>% 
    dplyr::ungroup() %>% 
    tidyr::pivot_wider(names_from = "name", values_from = "value") %>% 
    dplyr::select(dplyr::contains("quant_pos"))

  # add quant pos columns to dataframe with survey obs and PPD samples
  survey_estimate_ppd <- dplyr::bind_cols(
    survey_estimate_ppd, 
    survey_estimate_ppd_dist
  )

  # calculate position of oos obs within ordered PPD (i.e. estimate of hist)
  oos_within_ppd <- function(x, CI_range) {
    sum(
      x >= 0.5 * (1 - CI_range) * N &
        x <= N - 0.5 * (1 - CI_range) * N
    ) / length(x)
  }
  
  oos_within_ppd_percent <- oos_within_ppd(
    survey_estimate_ppd$quant_pos_final, CI_range
  )
  print(paste0(
    "Percentage of survey points which fall within posterior predictive",
    " distribution at ",
    # 95%
    CI_range * 100,
    "%",
    " CI: ",
    round(oos_within_ppd_percent * 100, 2),
    "%"
  ))

  # TODO: Add option to include additional summary statistics
  summary_stats <- list(
    "oos_observations_within_PPD_CI" = oos_within_ppd_percent
  )


  #### Calculate ELPD, CRPS & Error Stats ####

  # compute summary stats for comparison with other models, if specified
  if (compare_stats == TRUE) {
    
    # Actual OOS survey obs vs posterior predictive distribution for each
    actual <- survey_estimate_ppd$mean
    predictions <- dplyr::select(survey_estimate_ppd, dplyr::contains("samp_"))
    
    # calculate ELPD
    elpd <- loo::elpd(t(predictions))
    
    # calculate CRPS
    crps <- scoringutils::crps_sample(
      true_values = actual,
      predictions = as.matrix(predictions)
    )
    
    # calculate error stats (MAE, MSE & RMSE)
    
    # errors for each sample from PPD for each observation
    errors <- actual - as.matrix(predictions) # errors
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
