#' @title Produce Population weighted Aggregated Samples
#' @description Aggregate by area, year, age and type (weighted by population),
#' and convert to a percentage/probability.
#' @param .data \code{data.frame} including area populations, with
#' un-aggregated samples.
#' @param aggr_cols Columns to aggregate samples by, Default:
#' c("area_id", "area_name", "year", "age", "age_group", "model", "type")
#' @return \code{data.frame} with samples aggregated by \code{aggr_cols} and
#' weighted by population.
#' @importFrom dplyr %>%
#' @export
aggregate_sample <- function(.data,
                             aggr_cols = c(
                               "area_id", "area_name", "year",
                               "age", "age_group", "model", "type"
                             )) {

  # ensure aggregation columns are in the data
  aggr_cols <- aggr_cols[aggr_cols %in% names(.data)]

  # Multiplying by population to population weight
  .data <- .data %>%
    dplyr::mutate(dplyr::across(dplyr::contains("samp_"), ~ . * population))

  # summarise samples by aggr_cols
  .data <- data.table::setDT(.data)[,
    lapply(.SD, sum, na.rm = T),
    by = c(aggr_cols),
    .SDcols = c("population", paste0("samp_", c(1:100)))
  ]
  # divide by population to population weight
  .data <- as.data.frame(.data) %>%
    dplyr::mutate(dplyr::across(dplyr::contains("samp_"), ~ . / population))

  return(.data)
}
