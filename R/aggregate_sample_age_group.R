#' @title Produce Population weighted Aggregated Samples for Age Groups
#' @description Aggregate by area, year, age group (rather than discrete ages)
#' and type (weighted by population), and convert to a percentage/probability.
#' @param results_list list of \code{data.frame}s outputted by
#' \code{\link[threemc]{combine_areas}} with \code{join = FALSE}, including
#' area populations, with un-aggregated samples.
#' @param aggr_cols Columns to aggregate samples by, Default:
#' c("area_id", "area_name", "year", "model", "type")
#' @param age_groups Age groups to aggregate by, Default:
#' c("0-4",   "5-9",   "10-14", "15-19", "20-24", "25-29",
#' "30-34", "35-39", "40-44", "45-49", "50-54", "54-59",
#' "0+",    "10+",   "15+",   "15-24", "10-24", 15-29",
#' "10-29", "15-39", "10-39", "15-49", "10-49")
#' @param N Number of samples to summarise, Default: NULL
#' @return \code{data.frame} with samples aggregated by \code{aggr_cols} and
#' weighted by population.
#' @seealso
#'  \code{\link[threemc]{combine_areas}}
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
aggregate_sample_age_group <- function(results_list,
                                       aggr_cols = c("area_id", "area_name", 
                                                     "year", "model", "type"),
                                       age_groups = c(
                                         "0-4", "5-9", "10-14", "15-19", 
                                         "20-24", "25-29", "30-34", "35-39", 
                                         "40-44", "45-49", "50-54", "54-59",
                                         "0+", "10+", "15+", "15-24", "10-24",
                                         "15-29", "10-29", "15-39", 
                                         "10-39", "15-49", "10-49"
                                       ),
                                       N = 100
                                       ) {
  if (inherits(results_list, "data.frame")) {
    stop("requires list from combine_areas (set argument join = FALSE)")
  }

  #global bindings for data.table non-standard evaluation
  .SD <- NULL
  
  # Multiplying by population to population weight
  results_list <- lapply(results_list, function(x) {
    x %>%
      dplyr::mutate(
        dplyr::across(dplyr::contains("samp_"), ~ . * population)
      )
  })
  # aggregate sample for each age group
  results <- lapply(seq_along(age_groups), function(i) {
    # If upper limit use this split
    if (grepl("-", age_groups[i])) {
      age1 <- as.numeric(strsplit(age_groups[i], "-")[[1]][1])
      age2 <- as.numeric(strsplit(age_groups[i], "-")[[1]][2])
    }
    # If no upper limit use this split
    if (grepl("\\+", age_groups[i])) {
      age1 <- as.numeric(strsplit(age_groups[i], "\\+")[[1]][1])
      age2 <- Inf
    }
    results_list_loop <- lapply(results_list, function(x) {
      x <- x %>%
        # take results for age group i
        dplyr::filter(.data$age >= age1, .data$age <= age2) %>%
        dplyr::select(-.data$age)
      # Getting summarising samples
      x <- data.table::setDT(x)[,
        lapply(.SD, sum, na.rm = T),
        by = c(aggr_cols),
        .SDcols = c("population", paste0("samp_", c(1:N)))
      ]
      x <- x %>%
        # Adding age group
        dplyr::mutate(age_group = age_groups[i])
    })
    # Printing index
    print(age_groups[i])
    # return ages
    return(results_list_loop)
  })
  # join together
  results <- as.data.frame(data.table::rbindlist(
    lapply(results, data.table::rbindlist)
  ))

  # Multiplying by population to population weight
  # (don"t do this for "N performed", if present)
  results <- results %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::contains("samp_"), ~ ifelse(grepl("performed", type),
          .,
          . / population
        )
      )
    )

  return(results)
}
