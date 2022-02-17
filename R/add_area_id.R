#' @title 
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
#' @export
add_area_id <- function(
  df, 
  df_areas_wide, 
  par, 
  add_keep_cols = NULL) {
  
    # Getting area_id's
    area_lev_current_id = paste0("area_id", par$area_lev)
    # The level we want
    area_lev_select_id = paste0("area_id", par$area_lev_select)
    area_lev_select_name = paste0("area_name", par$area_lev_select)

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
                      distinct(),
                  by = c("area_id")) %>%
        # Select the right columns (account for case when we are at the lowest level)
        dplyr::select(
            "area_id" = ifelse(par$area_lev_select == par$area_lev,
                               "area_id",
                               area_lev_select_id),
            "area_name",
            dplyr::all_of(select_cols),
            dplyr::all_of(par$sample_cols)
        )
    return(df_area_id)
}
