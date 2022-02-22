#' @title Collect results for lower area hierarchies
#' @description Function to collect results for lower area hierarchies by
#' joining higher area hierarchies.
#' @param .data Results for highest area hierarchy, to be combined to
#' give results for lower/less granular area hierarchies.
#' @param areas_wide \code{data.frame} with shapefiles and area hierarchy.
#' @param area_lev Desired admin boundary level.
#' @param join Indicator to decide whether to join data for different
#' area hierarchies, or return them in list form.
#' @param add_keep_cols Additional columns to keep when summarising,
#' Default: NULL
#' @param ... Further arguments passed to \link[threemc]{add_area_id}.
#' @return \code{data.frame} or list (depending on the value of \code{join})
#' with results for all area levels less than or equal to \code{area_lev}.
#' @importFrom dplyr %>%
#' @importFrom data.table %like%
#' @importFrom rlang .data
#' @export

# hierarchies
combine_areas <- function(.data,
                          areas_wide,
                          area_lev,
                          join,
                          add_keep_cols = NULL,
                          ...) {

  # all area levels in the data (0 indexed)
  area_levs <- seq_len(area_lev) - 1

  # columns to keep
  add_keep_cols <- c(add_keep_cols, names(.data)[names(.data) %like% "samp_"])
  if (length(add_keep_cols) == 0) add_keep_cols <- NULL

  # collect results for lower area hierarchies by joining higher area
  # hierarchies (do so until you reach "area 0")
  if (area_levs[1] > -1) {
    results_list <- lapply(area_levs, function(x) {
      add_area_id(
        df = (.data %>% dplyr::select(-.data$area_name)),
        df_areas_wide = areas_wide,
        par = list(
          "area_lev" = area_lev,
          "area_lev_select" = x
        ),
        add_keep_cols = add_keep_cols
      )
    })
    # also add highest area hierarchy data
    results_list[[length(results_list) + 1]] <- .data
  } else {
    results_list <- .data
  }

  # return list or dataframe?
  if (join) {
    return(as.data.frame(data.table::rbindlist(
      results_list,
      use.names = T, ...
    )))
  } else {
    return(results_list)
  }
}
