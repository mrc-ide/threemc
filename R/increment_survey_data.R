#' @title Increase Area Level
#' @description Function to increase area levels.
#' @param survey_data Survey data, whose area level in \code{area_id} is not as 
#' desired.
#' @param areas_wide \code{data.frame} with shapefiles and area hierarchy.
#' @param par list with two entries: 
#' \itemize{
#'  \item{\code{area_lev}}{Current area level of \code{df}.}
#'  \item{\code{area_lev_select}}{Desired area level for \code{df}.}
#' }
#' @return Survey data, with desired area level in \code{area_id}.
#' @importFrom dplyr %>%
#' @export
# levels (should work for decreasing as well??)
increment_survey_area <- function(survey_data,
                                  areas_wide,
                                  par) {

    # store to later keep only these
    orig_names <- names(survey_data)

    # take only those areas in the area level you want increased:
    survey_data_area_lev <- survey_data %>%
        dplyr::mutate(area_level = as.numeric(substr(area_id, 5, 5))) %>%
        dplyr::filter(area_level == !!par$area_lev)
    survey_data <- survey_data %>%
        dplyr::anti_join(survey_data_area_lev)

    # change area level to desired level
    # (should expand this to use combine_areas from aggregations!)
    survey_data_area_lev <- add_area_id(
        df = (survey_data_area_lev),
        df_areas_wide = areas_wide,
        par = par,
        add_keep_cols = names(survey_data_area_lev)
    ) %>%
        dplyr::select(dplyr::all_of(orig_names))

    # join back with other surveys
    survey_data <- rbind(survey_data_area_lev, survey_data)
}
