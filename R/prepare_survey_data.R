#' @title Prepare Survey Data
#' @description Prepare survey data required to run the circumcision model. Can
#' also optionally apply \link[threemc]{normalise_weights_kish}, to
#'  normalise survey weights and apply Kish coefficients.
#' @param areas `sf` shapefiles for specific country/region.
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
#' @param strata.norm Stratification variables for normalising survey weights.
#' @param strata.kish Stratification variables for estimating and applying the
#' Kish coefficients.
#' @seealso
#'  \code{\link[threemc]{normalise_weights_kish}}
#' @return Survey data with required variables to run circumcision model.
#' @export
#'
#' @import dplyr
#' @import sf
#' @import rlang

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

  # change colnames to those in line with areas
  if ("geoloc_area_id" %in% names(survey_clusters)) {
      survey_clusters <- rename(survey_clusters, area_id = geoloc_area_id)
  }

  ## Bringing datasets together
  survey_circumcision <- survey_circumcision %>%
    ## Merging on individual information to  the circumcision dataset
    left_join(
      (survey_individuals %>%
        select(contains("id"), any_of(c("sex", "age", "indweight")))),
      by = c("survey_id", "individual_id")
    ) %>%
      ## Merging on cluster information to the circumcision dataset
    left_join(
      (survey_clusters %>%
       select(any_of(c("survey_id", "cluster_id")),
              "area_id" = "area_id")),
      by = c("survey_id", "cluster_id")
    ) %>%
    ## Remove those with missing circumcison status
    filter(
        !is.na(.data$circ_status),
        !is.na(.data$age),
        !is.na(.data$indweight)
    ) %>%
    ## Variables needed for analysis
    mutate(
      ## Survey year
      year = as.numeric(substr(.data$survey_id, 4, 7)),
      ## Year of Birth (estimated as no DOB filly yet)
      yob = .data$year - .data$age,
      ## If circumcision age > age of the individual set, reset circumcision age
      circ_age = ifelse(.data$circ_age > .data$age, NA, .data$circ_age)
    )

  ## Censoring if necessary ---------------------------------------------------

  ## Censoring at cens_year if assumed no circumcisions after a certain year
  if (!is.null(cens_age)) {
    survey_circumcision <- survey_circumcision %>%
      ## Censoring individuals from analysis at cens_age
      mutate(
        ## No circumcision after cens_age
        circ_status = ifelse(.data$circ_status == 1 &
          !is.na(.data$circ_age) &
          .data$circ_age > cens_age, 0.0, .data$circ_status),
        ## Resetting age at circumcision
        circ_age = ifelse(.data$circ_age > cens_age, NA,
                                .data$circ_age),
        ## Resetting age for everyone else
        age = ifelse(.data$age > cens_age, cens_age,
                           .data$age),
        ## Year of circ/censoring (estimated using the age as no date of circ)
        yoc = ifelse(!is.na(.data$circ_age), .data$yob + .data$circ_age,
                           .data$yob + .data$age)
      )
  }

  ## Censoring at cens_year if assumed no circumcisions after a certain year
  if (!is.null(cens_year)) {
    survey_circumcision <- survey_circumcision %>%
        ## Censoring at cens_year
        filter(.data$yob < cens_year) %>%
        ## Final variables for modelling
        mutate(
            ## Censoring circumcision status for those circumcised in cens_year,
            ## Assuming interval censored people were circumcised before cens_year
            circ_status = ifelse(.data$yoc >= cens_year &
                                     .data$circ_status == 1 & !is.na(.data$circ_age),
                                 0.0, .data$circ_status),
            ## circ censoring year / censor year in cens_year - 1 at cens_year - 1
            yoc = ifelse(.data$yoc == cens_year, cens_year - 1, .data$yoc)
        )
  }

  ## Setting desired level aggregation ----------------------------------------

  ## Getting the area level id to province
  for (i in seq_len(max(areas$area_level))) {
    survey_circumcision <- survey_circumcision %>%
      ## Merging on boundary information
      left_join(
        (areas %>%
          st_drop_geometry() %>%
          select(contains("area_id"), matches("area_level"))),
        by = "area_id"
      ) %>%
      ## Altering area
      mutate(
        area_id = if_else(.data$area_level == area_lev,
                          as.character(.data$area_id),
                          as.character(.data$parent_area_id))
        ) %>%
      select(-any_of(c("parent_area_id", "area_level")))
  }

  ## Final preparation of circumcision variables ------------------------------

  ## Preparing circumcision variables for the model
  survey_circumcision <- survey_circumcision %>%
    ## Merging on the region index
    left_join(
      (areas %>%
        st_drop_geometry() %>%
        select(any_of(c("area_id", "area_name", "space")))),
      by = "area_id"
    ) %>%
    mutate(
      ## Time interval for the individual
      time1 = .data$yob - start_year + 1,
      time2 = .data$yoc - start_year + 1,
      ## Event type
      event = ifelse(.data$circ_status == 1 & !is.na(.data$circ_age), 1,
        ifelse((.data$circ_status == 1 & is.na(.data$circ_age)), 2, 0)
      ),
      ## Circumcision age
      circ_age = .data$yoc - .data$yob,
      age = .data$circ_age + 1
    )

  ## Adding circumcision type to dataset
  survey_circumcision <- survey_circumcision %>%
    ## Type of circumcision
    mutate(
      circ_who = ifelse(.data$circ_who == "other",
        NA_character_,
        .data$circ_who
      ),
      circ_where = ifelse(.data$circ_where == "other",
        NA_character_,
        .data$circ_where
      ),
      type = case_when(
        .data$circ_who == "medical" | .data$circ_where == "medical" ~ "MMC",
        .data$circ_who == "traditional" | .data$circ_where == "traditional" ~ "TMC",
        TRUE ~ "Missing"
      )
    )

  ## Getting surveys without any type information
  tmp <- with(survey_circumcision, as.data.frame(table(survey_id, type))) %>%
    group_by(.data$survey_id) %>%
    ## calculate percentage and find surveys with all missing data
    mutate(Freq = .data$Freq / sum(.data$Freq)) %>%
    filter(.data$type == "Missing", .data$Freq == 1)

  ## return message detailing all surveys which are missing
  if (nrow(tmp) > 1) {
    for (i in seq_len(nrow(tmp))) {
      message(paste0(tmp$survey_id[i], " has all type == \"Missing\", and will
                     be removed from the data"))
    }
  }

  ## Removing surveys and individuals without any type information
  survey_circumcision <- survey_circumcision %>%
    filter(
      !(.data$survey_id %in% !!tmp$survey_id),
      !(.data$circ_status == 1 & .data$type == "Missing"),
      !is.na(.data$space)
    )

  ## normalise survey weights and apply Kish coefficients, if desired
  if (norm_kisk_weights) {
    survey_circumcision <- threemc::normalise_weights_kish(
      survey_circumcision,
      strata.norm = strata.norm,
      strata.kish = strata.kish
    )
  }

  ## Returning prepped circumcision datasets
  return(survey_circumcision)
}
