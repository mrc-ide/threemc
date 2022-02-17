#' @title Normalise Survey Weights and apply Kish Coefficients
#'
#' @description Normalise survey weights and apply Kish coefficients.
#'
#' @param survey_circumcision Information on male circumcision status from
#' surveys containing survey weights.
#' @param strata.norm Stratification variables for normalising survey weights,
#' Default: c("survey_id", "area_id")
#' @param strata.kish Stratification variables for estimating and applying the
#' Kish coefficients, Default: "survey_id"
#'
#' @return Survey data with normalised survey weights and required variables to
#' run circumcision model.
#' @export
#'
#' @importFrom dplyr %>%
#' @import rlang
normalise_weights_kish <- function(survey_circumcision,
                                   strata.norm = c("survey_id", "area_id"),
                                   strata.kish = c("survey_id")) {

  ## Preparing survey weights for the model
  survey_circumcision <- survey_circumcision %>%
    ## Standardising survey weights
    dplyr::group_by(dplyr::across(dplyr::all_of(strata.norm))) %>%
    dplyr::mutate(
      indweight_st = .data$indweight / mean(.data$indweight, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    ## Applying Kish coefficient to the survey weights
    dplyr::left_join(
      (survey_circumcision %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(strata.kish))) %>%
        dplyr::summarise(
          N = length(.data$survey_id),
          Neff = (sum(.data$indweight)^2) /
            sum(.data$indweight * .data$indweight),
          ratio = .data$N / .data$Neff,
          .groups = "drop"
        )),
      by = "survey_id"
    ) %>%
    dplyr::mutate(indweight_st = .data$indweight_st / .data$ratio)
}