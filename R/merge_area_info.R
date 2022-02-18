#' @title Merge Regional Informatoin on Dataset
#' @description Merge regional information on the dataset
#' (i.e. parent area info).
#' @param results \code{data.frame} you wish to merge shapefiles with.
#' @param areas \code{sf} shapefiles for specific country/region.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
merge_area_info <- function(results, areas) {

  # Merging regional information on the dataset (i.e. parent area info)
  results <- results %>%
    # Adding region information
    dplyr::left_join(
      (areas %>%
        dplyr::select(.data$area_id:.data$area_level)),
      by = c("area_id", "area_name")
    ) %>%
    dplyr::relocate(.data$area_level, .after = "area_name") %>%
    dplyr::left_join(
      (areas %>%
        dplyr::select(
          parent_area_id = .data$area_id,
          parent_area_name = .data$area_name
        )),
      by = "parent_area_id"
    ) %>%
    dplyr::relocate(dplyr::contains("area"))
}
