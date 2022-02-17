# collect results for lower area hierarchies by joining higher area
# hierarchies (should really allow inputs to add_keep_cols here!)
combine_areas <- function(
    .data,
    areas_wide,
    area_lev,
    join,
    add_keep_cols = NULL,
    ...
) {

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
                df = (.data %>% select(-area_name)),
                df_areas_wide = areas_wide,
                par = list("area_lev" = area_lev,
                           "area_lev_select" = x),
                add_keep_cols = add_keep_cols
            )}
        )
        # also add highest area hierarchy data
        results_list[[length(results_list) + 1]] <- .data
    } else {
        results_list <- .data
    }

    # return list or dataframe?
    if (join == TRUE) {
        return(as.data.frame(data.table::rbindlist(
            results_list, use.names = T, ...)
        ))
    } else {
        return(results_list)
    }
}