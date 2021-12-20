#' @title Prepare Survey Data
#' 
#' @description Prepare survey data required to run the circumcision model. Can
#' also optionally apply \link[threemc]{normalise_weights_kish}, to 
#'  normalise survey weights and apply Kish coefficients. 
#' 
#' @param area_hierarchy - Hierarchy and metadata of administrative boundaries.
#' @param area_boundaries - SF file with administrative boundaries (shapefiles).
#' @param survey_circumcision - Information on male circumcision status from 
#' surveys.
#' @param survey_individuals - Information on the individuals surveyed.
#' @param survey_clusters - Information on the survey clusters.
#' @param area_lev - Desired admin boundary level to perform the analysis on.
#' @param start_year - Year to begin the analysis on.
#' @param cens_year - Year to censor the circumcision data by (Sometimes some 
#' weirdness at the final survey year, e.g. v small number of MCs).
#' @param cens_age - Age to censor the circumcision data at (Default 59, 
#' i.e. no circumcisions after 59 years of age).
#' @param norm_kisk_weights - Set == TRUE to normalise survey weights and 
#' apply Kish coefficients.
#' @param strata.norm - See help file for \link[threemc]{normalise_weights_kish}.
#' @param strata.kish - See help file for \link[threemc]{normalise_weights_kish}.
#' 
#' @return Survey data with required variables to run circumcision model.
#' @export
#' 
#' @import dplyr
#' @import sf
#' @import threemc
prepare_survey_data <- function(areas,
                                survey_circumcision,
                                survey_individuals,
                                survey_clusters,
                                area_lev = 0,
                                start_year = 2006,
                                cens_year = NULL,
                                cens_age = 59,
                                norm_kisk_weights = F,
                                strata.norm = c("survey_id", "area_id"),
                                strata.kish = c("survey_id")) {
  
  ## Merging circumcision and individuals survey datasets ---------------------
  
  # Bringing datasets together
  survey_circumcision <- survey_circumcision %>%
    # Merging on individual information to  the circumcision dataset
    left_join(
      (survey_individuals %>%
        select(contains("id"), sex, age, indweight)),
      by = c("survey_id", "individual_id")) %>%
    # Merging on cluster information to the circumcision dataset
    left_join(
      (survey_clusters %>%
        select(survey_id, cluster_id, area_id = geoloc_area_id)),
      by = c("survey_id", "cluster_id")) %>%
    # Remove those with missing circumcison status
    filter(!is.na(circ_status), !is.na(age), !is.na(indweight)) %>%
    # Variables needed for analysis
    mutate(
      # Survey year
      year = as.numeric(substr(survey_id, 4, 7)),
      # Year of Birth (estimated as no DOB filly yet)
      yob = year - age,
      # If circumcision age > age of the individual set, reset circumcision age
      circ_age = ifelse(circ_age > age, NA, circ_age)
    )
  
  ## Censoring if necessary ---------------------------------------------------
  
  # Censoring at cens_year if assumed no circumcisions after a certain year
  if (!is.null(cens_age)) {
    survey_circumcision <- survey_circumcision %>%
      # Censoring individuals from analysis at cens_age
      mutate(
        # No circumcision after cens_age
        circ_status = ifelse(circ_status == 1 & !is.na(circ_age) & 
                               circ_age > cens_age, 0.0, circ_status),
        # Resetting age at circumcision
        circ_age = ifelse(circ_age > cens_age, NA, circ_age),
        # Resetting age for everyone else
        age = ifelse(age > cens_age, cens_age, age),
        # Year of circ/censoring (estimated using the age as no date of circ)
        yoc = ifelse(!is.na(circ_age), yob + circ_age, yob + age)
      )
  }
  
  # Censoring at cens_year if assumed no circumcisions after a certain year
  if (!is.null(cens_year)) {
    survey_circumcision <- survey_circumcision %>%
      # Censoring at cens_year
      filter(yob < cens_year) %>%
      # Final variables for modelling
      mutate(
        # Censoring circumcision status for those circumcised in cens_year,
        # Assuming the interval censored people were circumcised before cens_year
        circ_status = ifelse(yoc >= cens_year & circ_status == 1 & 
                               !is.na(circ_age), 0.0, circ_status),
        # circ censoring year (or censor year in cens_year - 1) at cens_year - 1
        yoc = ifelse(yoc == cens_year, cens_year - 1, yoc)
      )
  }
  
  ## Setting desired level aggregation ----------------------------------------
  
  # Getting the area level id to province
  for (i in seq_len(max(areas$area_level))) {
    survey_circumcision <- survey_circumcision %>%
      # Merging on boundary information
      left_join(
        (areas %>%
          st_drop_geometry() %>%
          select(contains("area_id"), area_level)),
        by = "area_id"
      ) %>%
      # Altering area
      mutate(area_id = if_else(area_level == area_lev,
                               as.character(area_id),
                               as.character(parent_area_id))) %>%
      select(-c(parent_area_id, area_level))
  }
  
  ## Final preparation of circumcision variables ------------------------------
  
  # Preparing circumcision variables for the model
  survey_circumcision <- survey_circumcision %>%
    # Merging on the region index
    left_join(
      areas %>%
        st_drop_geometry() %>%
        select(area_id, area_name, space),
      by = "area_id"
    ) %>%
    mutate(
      # Time interval for the individual
      time1 = yob - start_year + 1,
      time2 = yoc - start_year + 1,
      # Event type
      event = ifelse(circ_status == 1 & !is.na(circ_age), 1,
                     ifelse((circ_status == 1 & is.na(circ_age)), 2, 0)),
      # Circumcision age
      circ_age = yoc - yob,
      age = circ_age + 1)
  
  # Adding circumcision type to dataset
  survey_circumcision <- survey_circumcision %>%
    # Type of circumcision
    mutate(
      circ_who = ifelse(circ_who == "other", 
                        NA_character_, 
                        circ_who),
      circ_where = ifelse(circ_where == "other", 
                          NA_character_, 
                          circ_where),
      type = case_when(
        circ_who == "medical" | circ_where == "medical"         ~ "MMC",
        circ_who == "traditional" | circ_where == "traditional" ~ "TMC",
        TRUE                                                    ~ "Missing"
      )
    )
  
  # Getting surveys without any type information 
  tmp <- with(survey_circumcision, as.data.frame(table(survey_id, type))) %>% 
    group_by(survey_id) %>%
    # calculate percentage and find surveys with all missing data
    mutate(Freq = Freq / sum(Freq)) %>%
    filter(type == "Missing", Freq == 1)
  
  # return warning detailing all surveys which are missing
  if (nrow(tmp) > 1) {
    for (i in seq_len(nrow(tmp))) {
      warning(paste0(tmp$survey_id[i], " has all type == \"Missing\", and will 
                     be removed from the data"))
    }
  }

  # Removing surveys and individuals without any type information
  survey_circumcision <- survey_circumcision %>%
    filter(
      !(survey_id %in% !!tmp$survey_id),
      !(circ_status == 1 & type == "Missing"),
      !is.na(space)
    )
  
  # normalise survey weights and apply Kish coefficients, if desired
  if (norm_kisk_weights) {
    survey_circumcision <- normalise_weights_kish(survey_circumcision, 
                                                  strata.norm = strata.norm, 
                                                  strata.kish = strata.kish)
  }
  
  # Returning prepped circumcision datasets
  return(survey_circumcision)
}
