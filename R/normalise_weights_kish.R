#' @title Normalise Survey Weights and apply Kish Coefficients
#'
#' @description Normalise survey weights and apply Kish coefficients.
#'
#' @param survey_circumcision Information on male circumcision status from
#' surveys containing survey weights.
#' @param strata.norm Stratification variables for normalising survey weights.
#' @param strata.kish Stratification variables for estimating and applying the
#' Kish coefficients.
#'
#' @return Survey data with normalised survey weights and required variables to
#' run circumcision model.
#' @export
#'
#' @import dplyr
#' @import rlang
normalise_weights_kish <- function(survey_circumcision,
                                   strata.norm = c("survey_id", "area_id"),
                                   strata.kish = c("survey_id")) {

  ## Preparing survey weights for the model
  survey_circumcision <- survey_circumcision %>%
    ## Standardising survey weights
    group_by(across(all_of(strata.norm))) %>%
    mutate(
      indweight_st = .data$indweight / mean(.data$indweight, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    ## Applying Kish coefficient to the survey weights
    left_join(
      (survey_circumcision %>%
        group_by(across(all_of(strata.kish))) %>%
        summarise(
          N = length(.data$survey_id),
          Neff = (sum(.data$indweight)^2) /
            sum(.data$indweight * .data$indweight),
          ratio = .data$N / .data$Neff,
          .groups = "drop"
        )),
      by = "survey_id"
    ) %>%
    mutate(indweight_st = .data$indweight_st / .data$ratio)
}
