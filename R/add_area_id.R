#' change area ids from one hierarchy level to another
#' (works going "down" (i.e. less granular), but does vice versa work?)
add_area_id <- function(df, df_areas_wide, par, add_keep_cols = NULL) {

    # Getting area_id's
    area_lev_current_id = paste0("area_id", par$area_lev)
    # The level we want
    area_lev_select_id = paste0("area_id", par$area_lev_select)
    area_lev_select_name = paste0("area_name", par$area_lev_select)

    #' only select columns in our dataframe ("model" may be missing,
    #' and `age` and `age_group` are interchangable (can definitely improve this!!)
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
        left_join(df_areas_wide %>%
                      dplyr::select(
                          area_id = all_of(area_lev_current_id), # current level
                          all_of(area_lev_select_id), # desired level and
                          area_name = all_of(area_lev_select_name) # corresponding name
                      ) %>%
                      distinct(),
                  by = c("area_id")) %>%
        # Select the right columns (account for case when we are at the lowest level)
        dplyr::select(
            "area_id" = ifelse(par$area_lev_select == par$area_lev,
                               "area_id",
                               area_lev_select_id),
            "area_name",
            all_of(select_cols),
            all_of(par$sample_cols)
        )
    return(df_area_id)
}