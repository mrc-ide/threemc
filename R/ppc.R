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
threemc_oos_ppc <- function(
    fit,
    out,
    populations,
    survey_estimate,
    removed_years,
    type = "coverage",
    area_lev = 1,
    age_groups = c(
      # five-year age groups
      "0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
      "30-34", "35-39", "40-44", "45-49", "50-54", "54-59",
      # age groups with only minimum cut-off
      "0+", "10+", "15+",
      # other, wider age groups of interest
      "10-24", "15-24", "10-29", "15-29",
      "10-39", "15-39", "10-49", "15-49"
    ),
    CI_range = 0.95,
    N = 1000,
    compare_stats = TRUE
  ) {

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
    dplyr::filter(.data$year %chin% removed_years, .data$area_level == area_lev)


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


  #### Prepare Survey Points & Join with Samples ####

  survey_estimate_prep <- survey_estimate %>%
    # filter for OOS year(s) and modelled area_level
    dplyr::filter(
      .data$year %chin% removed_years,
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
    dplyr::summarise(mean = sum(.data$mean), .groups = "drop")

  # join with samples
  survey_estimate_ppd <- dplyr::left_join(
    survey_estimate_prep, out_types_agegroup
  )

  # check for NAs in first sample column
  stopifnot(!all(is.na(survey_estimate_ppd$samp_1)))

  
  #### Calculate Posterior Predictive Check for Prevalence Estimations ####

  # func calculating where in model sample distn empirical values are
  quant_pos_sum <- function(y, x) if (y < x) 0 else 1

  survey_estimate_ppd_dist <- survey_estimate_ppd %>%
    dplyr::relocate(.data$mean, .before = .data$samp_1) %>%
    dplyr::group_by(dplyr::across(.data$area_id:.data$area_level)) %>%
    dplyr::rowwise() %>%
    dplyr::summarise(
      quant_pos = sum(
        dplyr::across(dplyr::starts_with("samp_"), ~ quant_pos_sum(
          mean, .x
        ))
      ),
      .groups = "drop"
    )

  # calculate position of oos obs within ordered PPD (i.e. estimate of hist)
  oos_pos_ppd <- survey_estimate_ppd_dist$quant_pos
  # calculate percentage of oos obs within {CI_range}% CI of PPD
  oos_within_ppd_percent <- sum(
    oos_pos_ppd >= 0.5 * (1 - CI_range) * N &
      oos_pos_ppd <= N - 0.5 * (1 - CI_range) * N
  ) / length(oos_pos_ppd)
  
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
  
  
  #### Calculate ELPD and CRPS ####
  
  # compute summary stats for comparison with other models, if specified
  if (compare_stats == TRUE) {
    # calculate ELPD
    elpd <- loo::elpd(
      t(dplyr::select(survey_estimate_ppd, dplyr::contains("samp_")))
    )
    # calculate CRPS
    crps <- scoringutils::crps_sample(
      true_values = survey_estimate_ppd$mean, 
      predictions = as.matrix(
        dplyr::select(survey_estimate_ppd, dplyr::contains("samp_"))
      )
    )
    summary_stats <- c(summary_stats, list("elpd" = elpd, "crps" = crps))
  }

  # Return results
  return(list(
    "ppd"           = survey_estimate_ppd,
    "ppd_dist"      = survey_estimate_ppd_dist,
    "summary_stats" = summary_stats
  ))
}
