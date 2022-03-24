#' @title Create a List containing the Hierarchy of levels
#'
#' @description Create a list containing all the area dependencies and number
#' of for each area in the hierarchy
#'
#' @param areas `sf` shapefiles for specific country/region.
#' @param area_lev  PSNU area level for specific country.
#' @returns A list of length 2 containing:
#' \itemize{
#'  \item{"sub_region_list"}{A list of the specific sub-regions contained
#'  within each space (i.e. for each area_id) (including itself)}
#'  \item{"n_links_df"}{A dataframe with 2 columns, area_id and ndep, detailing
#'  the number of sub-regions contained within each area_id (also including
#'  itself)}
#' }
#' @rdname create_aggregate_structure
#' @export
create_aggregate_structure <- function(areas,
                                       area_lev) {
  # drop geometry and filter to specified area level
  areas <- sf::st_drop_geometry(areas) %>%
    dplyr::filter(.data$area_level <= area_lev)

  # Long to wide hierarchy (Need for new aggregation matrices)
  areas_wide <- spread_areas(areas)

  # iterate over all area_ids at specific area_lev (i.e. spaces)
  max_space <- areas %>%
    dplyr::filter(.data$area_level == area_lev) %>%
    dplyr::summarise(max(.data$space)) %>%
    dplyr::pull()
  area_id_seq <- seq(1, max_space, 1)
  sub_region_list <- lapply(area_id_seq, function(i) {
    # Getting areas lower in the hierarchy
    areas_wide %>%
      dplyr::filter(dplyr::if_any(dplyr::starts_with("space"), ~ . == i)) %>%
      dplyr::pull(paste0("space", area_lev))
  })

  n_sub_region_df <- areas %>%
    dplyr::distinct(.data$area_id) %>%
    dplyr::mutate(sp_dep = sapply(sub_region_list, length))

  # Returning list
  list(sub_region_list = sub_region_list, n_sub_region_df = n_sub_region_df)
}
