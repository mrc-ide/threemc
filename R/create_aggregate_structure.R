#' @title Create a List containing the Hierarchy of levels
#'
#' @description Create a list containing all the area dependencies and number 
#' of for each area in the hierarchy
#'
#' @param areas `sf` shapefiles for specific country/region.
#' @param area_lev  PSNU area level for specific country. 
#'
#' @rdname create_aggregate_structure
#' @export
create_aggregate_structure <- function(areas,
                                       area_lev) {
  # Long to wide hierarchy
  # Need this for the new aggregation matrices
  areas_wide <- areas %>%
    dplyr::filter(.data$area_level <= area_lev) %>%
    sf::st_drop_geometry() %>%
    spread_areas() 
  # Empty lists and data frame to store structure
  areas_agg1 <- list()
  areas_agg2 <- data.frame(
    area_id = subset(areas, area_level <= area_lev)$area_id, ndep = NA
  )
  # Loop for each area_id in the reference level
  for (i in seq_len((subset(areas, area_level == area_lev)$space))) {
    # Getting areas lower in the hierarchy
    test <- areas_wide %>%
      dplyr::filter(dplyr::if_any(dplyr::starts_with("space"), ~. == i)) %>%
      dplyr::pull(paste0("space", area_lev))
    # Adding list 
    areas_agg1 <- c(areas_agg1, list(test))
    areas_agg2$ndep[i] <- length(test)
  }
  # Returning list 
  return(list(areas_agg1 = areas_agg1, 
              areas_agg2 = areas_agg2))
}
