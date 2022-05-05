#' @title Aggregate survey points for each type and age group.
#' @description Aggregate survey points for each type and age group.
#' @param survey_circumcision Information on male circumcision status from
#' surveys.
#' @param areas \code{sf} shapefiles which include area hierarchies.
#' @param join Indicator to decide whether to join aggregated samples for
#' different age groups, Default: TRUE
#' @param types List of circumcision types to look at, Default:
#' list("Total" = c(unique(survey_circumcision$type)),
#'                 "Medical" = "MMC",
#'                 "Traditional" = "TMC")
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
#' @export
aggregate_sample_survey <- function(survey_circumcision,
                                    areas,
                                    join = TRUE,
                                    # types list
                                    types = list(
                                      "Total" = c(
                                        unique(survey_circumcision$type)
                                      ),
                                      "Medical" = "MMC",
                                      "Traditional" = "TMC"
                                    ),
                                    age_groups = c(
                                      "0-4", "5-9", "10-14", "15-19", "20-24",
                                      "25-29", "30-34", "35-39", "40-44",
                                      "45-49", "50-54", "54-59", "60-64", "65+",
                                      "0+", "10+", "15+", "15-24", "15-29",
                                      "15-39", "15-49", "10-29", "10-39",
                                      "10-49", "10-24"
                                    )) {

  # survey years
  survey_years <- unique(survey_circumcision$year)

  # loop through each type
  results_surv <- lapply(seq_along(types), function(i) {
    print(names(types)[i])

    # loop for each age group in the data (for each type)
    results_surv_type <- lapply(seq_along(age_groups), function(j) {
      print(age_groups[j])

      if (grepl("-", age_groups[j]) == TRUE) {
        age1 <- as.numeric(strsplit(age_groups[j], "-")[[1]][1])
        age2 <- as.numeric(strsplit(age_groups[j], "-")[[1]][2])
      }
      # If no upper limit use this split
      if (grepl("\\+", age_groups[j]) == TRUE) {
        age1 <- as.numeric(strsplit(age_groups[j], "+")[[1]][1])
        age2 <- Inf
      }
      # Getting proportions
      tmp <- survey_circumcision %>%
        dplyr::filter(.data$age >= age1, .data$age <= age2) %>%
        dplyr::group_by(.data$area_id, .data$year) %>%
        dplyr::summarise(
          Y_ind = length(.data$circ_status) *
            sum((.data$circ_status == 1 & .data$type %in% types[[i]]) *
              .data$indweight, na.rm = TRUE) /
            sum(.data$indweight, na.rm = TRUE),
          Y_obs = sum(.data$circ_status == 1 & .data$type %in% types[[i]]),
          N = length(.data$circ_status),
          p_ind = .data$Y_ind / .data$N,
          p_obs = .data$Y_obs / .data$N,
          .groups = "drop"
        ) %>%
        # adding age group
        dplyr::mutate(age_group = age_groups[i])
    })

    # append together results for each age group
    results_surv_type <- as.data.frame(
      data.table::rbindlist(results_surv_type, use.names = T)
    )

    # Adding to skeleton dataset and adding regional information
    results_surv_type <- expand.grid(
      area_id = sort(unique(results_surv_type$area_id)),
      year = survey_years,
      type = names(types)[i],
      age_group = age_groups
    ) %>%
      dplyr::left_join(results_surv_type) %>%
      # Adding region information
      dplyr::left_join(
        (areas %>%
          sf::st_drop_geometry() %>%
          dplyr::select(dplyr::contains("area"), -.data$area_level_label)),
        by = "area_id"
      ) %>%
      dplyr::left_join(
        (areas %>%
          sf::st_drop_geometry() %>%
          dplyr::select(
            parent_area_id = .data$area_id,
            parent_area_name = .data$area_name
          )),
        by = "parent_area_id"
      )
    return(results_surv_type)
  })

  if (join == TRUE) {
    results_surv <- as.data.frame(data.table::rbindlist(results_surv))
  }

  return(results_surv)
}
