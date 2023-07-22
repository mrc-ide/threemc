#' @title Posterior Predictive Distribution and checks on OOS survey
#' @description Aggregate specified `numeric` columns by population-weighted
#' age groups (rather than single year ages), split by specified categories. 
#' Using an alternative method to previously. 
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
#' @rdname threemc_ppc2
#' @export
threemc_ppc2 <- function(
    fit,
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
    # CI_range = c(0.5, 0.8, 0.95),
    N = 1000,
    seed = 123
  ) {
  
  # Take sample matrix rows that are kept in out from original skeleton data
  if ("n" %in% names(out)) {
    fit$sample <- lapply(fit$sample, function(x) x[out$n, ])
    out$n <- NULL
  } 
  
  # Extracting samples of probabilities from both models
  # Full == fit with program data, which we don't have!
  mc_prop <- 1.0 - fit$sample$surv_mc
  tmc_prop <- fit$sample$cum_inc_tmc
  mmc_prop <- fit$sample$cum_inc_mmc
  
  # check for type information
  type_info <- TRUE
  if (is.null(mmc_prop) && is.null(tmc_prop)) type_info <- FALSE
  
  # Assign row index to output frame -> identify rows in samples
  out_idx <- out %>%
    dplyr::select(area_id, year, age) %>%
    dplyr::mutate(idx = row_number())
  
  survey_obs <- survey_circumcision_test %>% 
    dplyr::rename(circ = circ_status)
  
  # impute type info with "Missing" (i.e. have non-NA)
  if (type_info == FALSE) {
    survey_obs$circ_who <- survey_obs$circ_where <- "Missing"
  }
  
  # Remove if missing circumcision status, age, weight or area_id
  survey_obs <- survey_obs %>%
    dplyr::filter(
      !is.na(circ),
      !is.na(age),
      !is.na(indweight),
      !is.na(area_id),
      !(circ == 1 & is.na(circ_where) & is.na(circ_who))
    )
  
  # Assign survey year
  survey_obs <- survey_obs %>%
    dplyr::mutate(
      year = as.numeric(substr(survey_id, 4, 7))
    ) %>% 
    # Subset to age < 60
    dplyr::filter(age < 60)
  
  # Assign circumcision type
  survey_obs <- survey_obs %>%
    dplyr::mutate(
      type = dplyr::case_when(
        circ_who == "medical" & circ_where == "medical"         ~ "MMC",
        circ_who == "medical" & circ_where == "traditional"     ~ "MMC",
        circ_who == "traditional" & circ_where == "medical"     ~ "MMC",
        circ_who == "traditional" & circ_where == "traditional" ~ "TMC",
        is.na(circ_who) & circ_where == "medical"               ~ "MMC",
        is.na(circ_who) & circ_where == "traditional"           ~ "TMC",
        circ_who == "medical" & is.na(circ_where)               ~ "MMC",
        circ_who == "traditional" & is.na(circ_where)           ~ "TMC",
        is.na(circ_who) & is.na(circ_where)                     ~ "Missing"
      )
    )
  
  # Roll-up area_id to district (level 2) (no need to do for KEN!)
  if (any(survey_obs$area_level > area_lev)) {
   for (i in seq_len(max(survey_obs$area_level))) {
     survey_obs <- survey_obs %>%
       dplyr::select(-matches("area_level")) %>% 
       dplyr::left_join(
         sf::st_drop_geometry(areas) %>%
           dplyr::select(area_id, area_level, parent_area_id),
         by = "area_id"
       ) %>%
       dplyr::mutate(
         area_id = ifelse(
           area_level == area_lev, 
           as.character(area_id), 
           as.character(parent_area_id)
         )
       ) %>%
       dplyr::select(-c(parent_area_id, area_level))
   }
  }
  
  # Assign output frame idx
  survey_obs <- survey_obs %>%
    dplyr::left_join(out_idx, by = c("area_id", "year", "age")) %>% 
    dplyr::filter(!is.na(idx)) # sometimes NAs here!
  
  #### Simulate observations from the predictive distribution ####
  # Get posterior probability of being circumcised (any type) for each obs
  survey_mc_prop <- mc_prop[survey_obs$idx, ]
  
  # Simulate circumcised / uncircumcised as binomial 0/1
  # Reshape as matrix
  survey_sim_mc <- rbinom(length(survey_mc_prop), 1, survey_mc_prop)
  survey_sim_mc <- matrix(survey_sim_mc, nrow = nrow(survey_mc_prop))
  
  # Get probability of being MMC given circumcised (MMC or TMC)
  if (type_info) {
    survey_mmctype_prop <- mmc_prop[survey_obs$idx, ] / 
      (tmc_prop[survey_obs$idx, ] + mmc_prop[survey_obs$idx, ])
  
    # Simulate whether MMC _if_ they were circumcised
    # NOTE: This is a bit inefficient because it samples for _all_ respondents
    #       rather than only those who are circumcised.
    survey_sim_mmctype <- rbinom(
      length(survey_mmctype_prop), 1, survey_mmctype_prop
    )
    survey_sim_mmctype <- matrix(
      survey_sim_mmctype, nrow = nrow(survey_mmctype_prop)
    )
    
    # Construct matrices for simulated MMC and TMC
    survey_sim_mmc <- survey_sim_mc == 1 & survey_sim_mmctype == 1
    survey_sim_mmc[] <- as.integer(survey_sim_mmc)
    survey_sim_tmc <- survey_sim_mc == 1 & survey_sim_mmctype == 0
    survey_sim_tmc[] <- as.integer(survey_sim_tmc)
  }
  
  # Construct data frame with all outputs for observed and simulated
  # MC, MMC, and TMC
  survey_sim_obs <- survey_obs %>%
    transmute(
      survey_id,
      individual_id,
      area_id,
      sex,
      age,
      year,
      indweight,
      idx,
      mc_obs = circ,
      mmc_obs = as.integer(circ == 1 & type == "MMC"),
      tmc_obs = as.integer(circ == 1 & type == "TMC"),
    )
  
  # Altering column names
  # TODO: Functionalise!
  colnames(survey_sim_mc) <- sprintf(
    "mc_sim%04d", seq_len(ncol(survey_sim_mc))
  )
  if (type_info) {
    colnames(survey_sim_mmc) <- sprintf(
      "mmc_sim%04d", seq_len(ncol(survey_sim_mmc))
    )
    colnames(survey_sim_tmc) <- sprintf(
      "tmc_sim%04d", seq_len(ncol(survey_sim_tmc))
    )
  }
  
  # Appending together with the observations
  survey_sim <- dplyr::bind_cols(
    survey_sim_obs, survey_sim_mc
  )
  if (type_info) {
    survey_sim <- dplyr::bind_cols(
     survey_sim, survey_sim_mmc, survey_sim_tmc
    )
  } 
  
  # Pivot longer with column for circumcision type
  survey_sim <- survey_sim %>%
    tidyr::pivot_longer(
      cols = c(starts_with("mc_"), starts_with("mmc_"), starts_with("tmc_")),
      names_pattern = "(.*)_(.*)",
      names_to = c("type", ".value")
    )
  
  if (type_info == FALSE) {
    survey_sim <- survey_sim %>% 
      filter(!type %in% c("mmc", "tmc"))
  }
  
  
  #### Summarise samples, get outputs ####
  # Helper function for calculating coverage
  calc_ppd_coverage <- function(interval, obs) {
    qlower <- quantile(
      dplyr::c_across(starts_with("sim")), 0.5 - interval / 2, names = FALSE
    )
    qupper <- quantile(
      dplyr::c_across(starts_with("sim")), 0.5 + interval / 2, names = FALSE
    )
    between(obs, qlower, qupper)
  }
  
  # Calculate coverage by survey x district x age 15-49 
  mccov_district_15to59 <- survey_sim %>%
    dplyr::filter(age %in% 15:59) %>%
    dplyr::summarise(
      dplyr::across(
        c(obs, starts_with("sim")), ~ stats::weighted.mean(.x, indweight)
      ),
      .by = c(type, survey_id, area_id)
    ) 
  
  ppd_district_15to59 <- mccov_district_15to59 %>%
    dplyr::filter(!is.na(obs)) %>% 
    dplyr::rowwise(type, survey_id, area_id, obs) %>%
    dplyr::summarise(
      ppd_mean = mean(dplyr::c_across(starts_with("sim"))),
      crps = ifelse(
        is.na(obs),
        NA,
        scoringutils::crps_sample(
          true_values = obs, 
          predictions = matrix(dplyr::c_across(starts_with("sim")), nrow = 1)
        )
      ),
      mae = mean(abs(dplyr::c_across(starts_with("sim")) - obs)),
      rmse = sqrt(mean((dplyr::c_across(starts_with("sim")) - obs) ^ 2)),
      CI.0.5 = calc_ppd_coverage(0.5, obs),
      CI.0.8 = calc_ppd_coverage(0.8, obs),
      CI.0.95 = calc_ppd_coverage(0.95, obs),
      mean = obs, # keep obs, want to use later
      .groups = "drop"
    )
  
  # Calculate coverage by survey x district x age 15-49 
  mccov_district_5year <- survey_sim %>%
    dplyr::mutate(
      age_group = cut(
        age, 
        0:12 * 5, 
        sprintf("%02d-%02d", 0:11 * 5, 0:11 * 5 + 4), 
        right = FALSE
      )
    ) %>%
    # dplyr::mutate(iso3 = cntry) %>% 
    dplyr::summarise(
      dplyr::across(
        c(obs, starts_with("sim")), ~ stats::weighted.mean(.x, indweight)
      ),
      .by = c(type, survey_id, area_id, age_group)
    ) 
  
  ppd_district_5year <- mccov_district_5year %>%
    dplyr::filter(!is.na(obs)) %>% 
    dplyr::rowwise(
      type, survey_id, area_id, age_group, obs
    ) %>%
    dplyr::summarise(
      ppd_mean = mean(dplyr::c_across(starts_with("sim"))),
      crps = ifelse(
        is.na(obs),
        NA,
        scoringutils::crps_sample(
          true_values = obs, 
          predictions = matrix(dplyr::c_across(starts_with("sim")), nrow = 1)
        )
      ),
      mae = mean(abs(dplyr::c_across(starts_with("sim")) - obs)),
      rmse = sqrt(mean((dplyr::c_across(starts_with("sim")) - obs) ^ 2)),
      CI.0.5 = calc_ppd_coverage(0.5, obs),
      CI.0.8 = calc_ppd_coverage(0.8, obs),
      CI.0.95 = calc_ppd_coverage(0.95, obs),
      mean = obs,
      .groups = "drop"
    )
  
  # PPD coverage by district x 15-49 years
  # ppd_district_15to59 %>%
  #   dplyr::summarise(
  #     dplyr::across(crps, sum),
  #     dplyr::across(c(mae, rmse, starts_with("CI")),  ~ mean(.x, na.rm = TRUE)),
  #     .by = c(type)) %>%
  #   identity()
  # "~/imperial_repos/threemc-orderly/ken_paed_cutoff_no_time_tmc_ppc_age.csv"
  
  # PPD coverage by district x 15-49 years
  ppd_district_5year <- ppd_district_5year %>%
    dplyr::summarise(
      dplyr::across(crps, sum),
      dplyr::across(
        c(mae, rmse, starts_with("CI"), mean, ppd_mean), 
        ~ mean(.x, na.rm = TRUE)
      ),
      .by = c(area_id, type, survey_id, age_group)
    )
  
  # have output matching original threemc_ppc
  ppd_district_5year <- ppd_district_5year %>% 
    dplyr::mutate(
      year = as.numeric(substr(survey_id, 4, 7)),
      type = paste(toupper(type), "coverage")  
    ) %>% 
    dplyr::select(-c(survey_id)) %>% 
    dplyr::relocate(
      area_id, 
      year,
      age_group, 
      type, 
      mean,
      ppd_mean, 
      ppd_0.500 = CI.0.5, 
      ppd_0.800 = CI.0.8,
      ppd_0.950 = CI.0.95,
      crps, 
      mae, 
      rmse
    )
  
  # calculate elpd
  elpd <- loo::elpd(t(ppd_district_5year$ppd_mean))
  ppd_district_5year <- ppd_district_5year %>% 
    mutate(elpd = elpd$pointwise[, 1]) %>% 
    relocate(elpd, .before = "crps")
  
  # summarise within the function for each type
  ppd_district_5year_sum <- ppd_district_5year %>% 
  summarise(
    dplyr::across(dplyr::all_of(c("elpd", "crps")), sum),
    dplyr::across(
      c(mae, rmse, starts_with("ppd_0")),  ~ mean(.x, na.rm = TRUE)
    ),
    .by = c(type)
  )
  
  # print summary
  print(ppd_district_5year_sum) 
  
  return(list(
    "ppc_df"         = ppd_district_5year,
    "ppc_summary_df" = ppd_district_5year_sum
  ))
}
