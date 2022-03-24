#' @title Subset and Prepare Shell Dataset to Only One Admin Boundary 
#'
#' @description  Subset the shell dataset to a specific area of interest and 
#' resets the space counter. The output dataset can be used to set up model 
#' components on the specified administrative boundaries. 
#'
#' @param dat Shell dataset used for modelling
#' @param area_lev  PSNU area level for specific country. Defaults to the
#' Defaults to the maximum area level found in `dat` if not supplied.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
create_shell_dataset_area <- function(dat, area_lev = NULL) {
  
  if (is.null(area_lev)) {
    message(
      "area_lev arg missing, taken as maximum area level in shell dataset"
    )
    area_lev <- max(dat$area_level, na.rm = TRUE)
  }
  
  # Only doing the matrices on the specified aggregation
  dat <- dat %>%
    dplyr::filter(.data$area_level == area_lev) %>%
    # Resetting counter on space
    dplyr::mutate(space = .data$space - min(.data$space) + 1)
  
  ## Returning matrix
  return(dat)
}
