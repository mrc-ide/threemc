#' @title Function to read in Circumcision Data
#'
#' @description Function to read in circumcision data to fit model. Handles
#' csv with \code{\link[data.table]{fread}} (but outputs data as a
#' `data.frame`), and geographical data with code{\link[sf]{read_sf}} (for which
#' it also adds unique identifiers for each `area_level`).
#'
#' @param path Path to data.
#' @param filters Optional named vector, whose values dictate the values
#' filtered for in the corresponding column names. Only supports filtering for
#' one value for each column. default: NULL
#' @param ... Further arguments passed to or from other methods.
#'
#' @seealso
#'  \code{\link[data.table]{fread}}
#'  \code{\link[sf]{read_sf}}
#' @return relevant data set, filtered as desired.
#' @export
#'
#' @import dplyr
#' @import sf
#' @import rlang
read_circ_data <- function(path, filters = NULL, selected = NULL, ...) {

  ## maybe add a warning for missing "circ" columns for surveys?? And add
  ## in NAs in this situation (look at KEN for this)

  ## read in data, depending on file type
  cond <- tools::file_ext(path)%in% c("geojson", "shp", "shx")
  if (cond == T) {
    .data <- read_sf(path, ...)
  } else {
      # selection prior to loading is allowed by fread
      .data <- as.data.frame(data.table::fread(path, select = c(selected), ...))
  }

  ## if desired, recursively filter data with provided `filters` vector
  if (!is.null(filters)) {
    cols <- names(filters)
    vals <- as.vector(filters[seq_along(filters)])
    for (i in seq_along(filters)) {
      if (!cols[i] %in% names(.data)) next
      ## change col i to symbol (if present), evaluate corresponding filter
      .data <- filter(.data, !!rlang::sym(cols[i]) == vals[i])
    }
  }

  # Select specific columns, if desired (and present) (no need to do for fread)
  if (!is.null(selected) & cond == T) {
      .data <- .data %>%
          select(all_of(selected[selected %in% names(.data)]))
  }

  ## for areas, add unique identifier within Admin code and merge to boundaries
  if (inherits(.data, "sf")) {
    .data <- .data %>%
      group_by(.data$area_level) %>%
      mutate(space = row_number()) %>%
      ungroup()
  }
  return(.data)
}
