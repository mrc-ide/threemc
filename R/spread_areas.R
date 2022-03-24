#' Spread area hierarchy to wide format
#'
#' @param areas area hierarchy data.frame
#' @param min_level integer specifying the minimum level wanted
#' @param max_level integer specifying the maximum level wanted
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
spread_areas <- function(areas, 
                         min_level = min(areas$area_level), 
                         max_level = max(areas$area_level)) {
  
  if (inherits(areas, "sf")) {
    boundaries <- areas %>% dplyr::select(.data$area_id)
    areas <- sf::st_drop_geometry(areas)
  }
  else {
    boundaries <- NULL
  }
  stopifnot(min_level >= min(areas$area_level))
  stopifnot(max_level <= max(areas$area_level))
  areas_wide <- areas %>% 
    dplyr::filter(.data$area_level == min_level) %>% 
    dplyr::select(
      # What's happening here?
      `:=`(!!paste0("area_id", min_level), .data$area_id), 
      `:=`(!!paste0("area_name", min_level), .data$area_name), 
      `:=`(!!paste0("space", min_level), .data$space)
    )
  
  for (level in (min_level + 1):max_level) {
    areas_wide <- areas_wide %>% 
      dplyr::left_join(
        areas %>% 
          dplyr::filter(.data$area_level == level) %>% 
          dplyr::select(
            `:=`(!!paste0("area_id", level), .data$area_id),
            `:=`(!!paste0("area_name", level), .data$area_name), 
            .data$parent_area_id,
            `:=`(!!paste0("space", level), .data$space), 
          ), 
        by = stats::setNames(
          c("parent_area_id"), c(paste0("area_id", level - 1L))
        )
      )
  }
  areas_wide$area_id <- areas_wide[[paste0("area_id", max_level)]]
  if (!is.null(boundaries)) {
    areas_wide <- sf::st_as_sf(
      dplyr::left_join(areas_wide, boundaries, by = "area_id")
    )
  }
  return(areas_wide)
}
