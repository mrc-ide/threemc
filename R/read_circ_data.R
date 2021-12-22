#' @title Function to read in Circumcision Data
#' 
#' @description Function to read in circumcision data to fit model. Handles 
#' csv with \link[data.table]{fread} (but outputs data as a `data.frame`), and 
#' geographical data with \link[sf]{read_sf} (for which it also adds unique 
#' identifiers for each `area_level`).
#' 
#' @param path Path to data.
#' @param filters Optional named vector, whose values dictate the values 
#' filtered for in the corresponding column names. Only supports filtering for 
#' one value for each column. default: NULL
#' 
#' @seealso 
#'  \code{\link[data.table]{fread}} 
#'  \code{\link[sf]{read_sf}} 
#' @return relevant data set, filtered as desired. 
#' @export
#'
#' @import dplyr
#' @importFrom data.table fread
#' @importFrom sf read_sf
#' @importFrom rlang sym
  
# maybe add a warning for missing "circ" columns for surveys?? And add them in
# NAs in this situation (look at KEN for this)
read_circ_data <- function(path, filters = NULL) {
  
  # read in data, depending on file type
  if (grepl(".geojson", path)) {
    .data <- read_sf(path)
  } else .data <- as.data.frame(fread(path))
  
  # if desired, recursively filter data with provided `filters` vector
  if (!is.null(filters)) {
    cols <- names(filters)
    vals <- as.vector(filters[seq_along(filters)])
    for (i in seq_along(filters)) {
      if (!cols[i] %in% names(.data)) next
      # change col i to symbol (if present), evaluate corresponding filter
      .data <- filter(.data, !!sym(cols[i]) == vals[i])
    }
  }
  # for areas, add unique identifier within Admin code and merge to boundaries
  if ("sf" %in% class(.data)) {
    .data <- .data %>%
      group_by(area_level) %>%
      mutate(space = row_number()) %>%
      ungroup()
  }
  return(.data)
}
