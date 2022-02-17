#' @title Format circumcision survey data with age data reformatted. 
#' @description Format circumcision survey data with age data reformatted.  
#' @param areas_wide \code{data.frame} with shapefiles and area hierarchy.
#' @param survey_circumcision - Information on male circumcision status from
#' surveys.
#' @param survey_individuals - Information on the individuals surveyed.
#' @param survey_clusters - Information on the survey clusters.
#' @return Survey data, with age data reformatted, at highest/most granular 
#' area level.
#' @importFrom dplyr %>%
#' @export
# 
prepare_survey_aggregation <- function(
  areas_wide,
  survey_circumcision,
  survey_individuals,
  survey_clusters
) {
  
  survey_circumcision <- survey_circumcision %>%
    # allocate random number from 0-5 for 95 entries
    dplyr::mutate(
      circ_age = dplyr::case_when(
        circ_age == 95 ~ sample(seq(0, 5), n(), replace = TRUE),
        TRUE               ~ circ_age),
      # Correct status for those who have information on
      # circumcision but are missing circumcision status
      dplyr::circ_status = ifelse(is.na(circ_status) & !is.na(circ_age),
                           1,
                           ifelse(is.na(circ_status) & !is.na(circ_age),
                                  1,
                                  ifelse(is.na(circ_status) &
                                           !is.na(circ_where),
                                         1,
                                         circ_status)))
    ) %>%
    # merging individual weights to:
    # Merging on individual information to  the circumcision dataset
    dplyr::left_join(
      (survey_individuals %>%
         dplyr::select(dplyr::contains("id"), "sex", "age", "indweight"))
    ) %>%
    # Merging on cluster information to  the circumcision dataset
    left_join(
      (dplyr::survey_clusters %>%
         dplyr::select(dplyr::contains("id"), - survey_region_id) %>%
         dplyr::distinct())
    ) %>%
    # Remove those with missing circumcision status
    dplyr::filter(
      !is.na(circ_status),
      !(is.na(age) & is.na(circ_age)),
      !is.na(geoloc_area_id)
    ) %>%
    # Adding age group and type of circumcision
    dplyr::mutate(
      age_group = as.numeric(cut(age,
                                 breaks = c(seq(0, 65, by = 5), Inf),
                                 labels = 1:14,
                                 right = FALSE,
                                 include.lowest = TRUE)),
      year = as.numeric(substr(survey_id, 4, 7)),
      # setting survey types
      type = dplyr::case_when(
        circ_who ==  "Healthcare worker"     ~ "MMC",
        tolower(circ_who) ==  "medical"      ~ "MMC",
        tolower(circ_where) == "medical"     ~ "MMC",
        circ_who == "Traditional practioner" ~ "TMC",
        tolower(circ_who) == "traditional"   ~ "TMC",
        tolower(circ_where) == "traditional" ~ "TMC",
        circ_status != 0                     ~ "Missing",
        TRUE                                 ~ NA_character_
      )
    ) %>%
    # Altering column names
    dplyr::rename(area_id = geoloc_area_id)
  
  # add area id's:
  # find last area column (i.e. highest area level present)
  last_area_id <- dplyr::last(
    names(areas_wide)[names(areas_wide) %like% "area_id" & 
                        nchar(names(areas_wide)) > 7]
  )
  
  # add all area data up to last area hierarchy, selecting only required
  survey_circumcision <- survey_circumcision %>%
    dplyr::left_join(areas_wide, by = "area_id") %>%
    # append highest level areas to area_id col
    dplyr::mutate(area_id = eval(parse(text = last_area_id))) %>%
    # Remove other area columns but keep area level of interest
    dplyr::select(
      -c(names(areas_wide)[!names(areas_wide) %in% "area_id"],
         age_group) # not the same age_group as we usually use
    )
}
