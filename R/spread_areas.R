#' Spread area hierarchy to wide format
#'
#' @param areas area hierarchy data.frame
#' @param min_level integer specifying the minimum level wanted
#' @param max_level integer specifying the maximum level wanted
#' @param space whether to include "space" columns. Excluding these returns the 
#' same object as \code{naomi::spread_areas}, Default: TRUE
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
spread_areas <- function(areas,
                         min_level = min(areas$area_level),
                         max_level = max(areas$area_level),
                         space = TRUE) {
  if (inherits(areas, "sf")) {
    boundaries <- areas %>% dplyr::select(.data$area_id)
    areas <- sf::st_drop_geometry(areas)
  } else {
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
  
  # create safe sequence from min_level + 1 to max_level to loop over
  level_seq <- seq_len(max_level)
  if (length(level_seq) == 0) {
    level_seq <- 0
  } else {
    level_seq <- level_seq[level_seq >= (min_level + 1)]
  }
  
  if (all(level_seq != 0)) {
    for (level in level_seq) {
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
  }
  
  areas_wide$area_id <- areas_wide[[paste0("area_id", max_level)]]
  if (!is.null(boundaries)) {
    areas_wide <- sf::st_as_sf(
      dplyr::left_join(areas_wide, boundaries, by = "area_id")
    )
  }
  
  # removing "space" columns returns same object as naomi::spread_areas
  if (space == FALSE) {
    areas_wide <- areas_wide %>% 
      dplyr::select(-dplyr::contains("space"))
  }
  
  return(areas_wide)
}
