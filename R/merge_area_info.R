#' @title Merge Regional Informatoin on Dataset
#' @description Merge regional information on the dataset 
#' (i.e. parent area info).
#' @param results \code{data.frame} you wish to merge shapefiles with.
#' @param areas \code{sf} shapefiles for specific country/region.
#' @importFrom dplyr %>%
#' @export
merge_area_info <- function(results, areas) {

    # Merging regional information on the dataset (i.e. parent area info)
    results <- results %>%
        # Adding region information
        dplyr::left_join(
            (areas %>%
                 dplyr::select(area_id:area_level)),
            by = c("area_id", "area_name")
        ) %>%
        dplyr::relocate("area_level", .after = "area_name") %>%
        dplyr::left_join(
            (areas %>%
                 dplyr::select(parent_area_id = area_id,
                        parent_area_name = area_name)),
            by = "parent_area_id"
        ) %>%
        dplyr::relocate(dplyr::contains("area"))
}