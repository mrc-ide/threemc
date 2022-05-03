#### %||% ####

`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}

#### add_area_id #### 

#' @title Change \code{area_id} from one hierarchy level to another
#' @description Function to change \code{area_id} from one hierarchy level to
#' another.
#' @param df Dataframe with \code{area_id} column.
#' @param df_areas_wide \code{sf} \code{dataframe} with shapefiles and area
#' hierarchy.
#' @param par list with two entries:
#' \itemize{
#'  \item{\code{area_lev}}{Current area level of \code{df}.}
#'  \item{\code{area_lev_select}}{Desired area level for \code{df}.}
#' }
#' @param add_keep_cols Additional columns to keep when summarising,
#' Default: NULL
#' @return \code{df} with `area_id` changed to `area_lev_select`.
#' @importFrom dplyr %>%
#' @rdname add_area_id
#' @keywords internal
add_area_id <- function(df,
                        df_areas_wide,
                        par,
                        add_keep_cols = NULL) {

  # Getting area_id's
  area_lev_current_id <- paste0("area_id", par$area_lev)
  # The level we want
  area_lev_select_id <- paste0("area_id", par$area_lev_select)
  area_lev_select_name <- paste0("area_name", par$area_lev_select)

  # only select columns in our dataframe ("model" may be missing,
  # and `age` and `age_group` are interchangable)
  select_cols <- c("year", "age", "age_group", "population", "type", "model")
  select_cols <- select_cols[select_cols %in% names(df)]
  # additional columns to keep, if supplied
  if (!is.null(add_keep_cols)) {
    select_cols <- unique(c(select_cols, add_keep_cols))
  }
  # remove columns which interfere with select below
  select_cols <- select_cols[!select_cols %in% c("area_id", "area_name")]

  df_area_id <- df %>%
    # join in area names for chosen area_id
    dplyr::left_join(df_areas_wide %>%
      dplyr::select(
        # current level
        area_id = dplyr::all_of(area_lev_current_id),
        # desired level and
        dplyr::all_of(area_lev_select_id),
        # corresponding name
        area_name = dplyr::all_of(area_lev_select_name)
      ) %>%
      dplyr::distinct(),
    by = c("area_id")
    ) %>%
    # Select the right columns (account for when we are at the lowest level)
    dplyr::select(
      "area_id" = ifelse(par$area_lev_select == par$area_lev,
        "area_id",
        area_lev_select_id
      ),
      "area_name",
      dplyr::all_of(select_cols),
      dplyr::all_of(par$sample_cols)
    )
  return(df_area_id)
}


#### combine_areas ####

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
#' @rdname combine_areas
#' @keywords internal
combine_areas <- function(.data,
                          areas_wide,
                          area_lev,
                          join,
                          add_keep_cols = NULL,
                          ...) {

  # all area levels in the data (0 indexed)
  area_levs <- seq_len(area_lev) - 1
  if (length(area_levs) == 0) area_levs <- -1

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
