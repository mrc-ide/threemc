#' @title Calculate number of people circumcised
#' @description Calculate number of people circumcised (as well as unmet need).
#' @param results Results with samples for number of circumcisions performed
#' in each region.
#' @importFrom dplyr %>%
#'@importFrom rlang .data
#' @export
n_circumcised <- function(results) {
  
  # Getting number of circumcised men
  n_circ <- split(results, results$type)

  # get circumcised population by type
  n_circ_type <- lapply(n_circ, function(x) {
    x %>%
      dplyr::mutate(
        dplyr::across(dplyr::contains("samp_"), ~ . * population),
        type = paste0(
          "Number circumcised (",
          stringr::str_remove(.data$type, " coverage"),
          ")"
        )
      )
  })
  # also calculate unmet need
  n_circ_type[[length(n_circ_type) + 1]] <- results %>%
    dplyr::filter(.data$type == "MC coverage") %>%
    dplyr::mutate(
      dplyr::across(dplyr::contains("samp_"), ~ population * (1 - .)),
      type = "Unmet need"
    )

  # Append together
  results_n <- as.data.frame(data.table::rbindlist(n_circ_type, use.names = T))
}
