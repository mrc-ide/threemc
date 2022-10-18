#' @title Prepare Survey Data
#' @description Prepare survey data required to run the circumcision model. Can
#' also optionally apply \link[threemc]{normalise_weights_kish}, to
#'  normalise survey weights and apply Kish coefficients.
#' @param areas \code{sf} shapefiles for specific country/region.
#' @param survey_circumcision - Information on male circumcision status from
#' surveys. If this is a list or contains more than one country, the function
#' is performed for each country present, returning a list.
#' @param survey_individuals - Information on the individuals surveyed.
#' @param survey_clusters - Information on the survey clusters.
#' @param area_lev - Desired admin boundary level to perform the analysis on.
#' @param start_year - Year to begin the analysis on, Default: 2006
#' @param cens_year - Year to censor the circumcision data by (Sometimes some
#' weirdness at the final survey year, e.g. v small number of MCs),
#' Default: NULL
#' @param cens_age - Age to censor the circumcision data at, Default: 59
#' @param rm_missing_type - Indicator to decide whether you would like to keep
#' surveys where there is no MMC/TMC disinction. These surveys may still be
#' useful for determining MC levels, Default: FALSE
#' @param norm_kisk_weights - Indicator to decide whether to normalise survey
#' weights and apply Kish coefficients, Default: TRUE
#' @param strata.norm Stratification variables for normalising survey weights,
#' Default: c("survey_id", "area_id")
#' @param strata.kish Stratification variables for estimating and applying the
#' Kish coefficients, Default: "survey_id"
#' @seealso
#'  \code{\link[threemc]{normalise_weights_kish}}
#' @return Survey data with required variables to run circumcision model.
#' @export
#'
#' @importFrom rlang .data
#' @importFrom dplyr %>%
prepare_survey_data <- function(areas,
                                survey_circumcision,
                                survey_individuals = NULL,
                                survey_clusters = NULL,
                                area_lev,
                                start_year = 2006,
                                cens_year = NULL,
                                cens_age = 59,
                                rm_missing_type = FALSE,
                                norm_kisk_weights = TRUE,
                                strata.norm = c("survey_id", "area_id"),
                                strata.kish = c("survey_id")) {
  
  
  # function to parse country specific psnu area levels if desired
  select_area_lev <- function(area_lev, cntry, is_add_data_present) {
    
    # if area_lev is provided (i.e as a value, rather than a df), just return 
    if (!inherits(area_lev, "data.frame")) {
      return(area_lev)
    } 
    # if area_lev is a df (i.e. threemc::datapack_psnu_area_level) ...
    
    # find area_level for specified country from area_lev df
    area_lev <- area_lev %>%
      dplyr::filter(.data$iso3 == cntry) %>%
      dplyr::pull(.data$psnu_area_level)
    
    # if area_lev for specified cntry is present in df, return 
    if (length(area_lev) == 1) return(area_lev)
    
    # if area_level for cntry is missing, assume most common lev in surveys
    if (length(area_lev) == 0 && is_add_data_present) {
      # If survey_clusters is specified, find most common area level there
      area_lev_df <- survey_clusters
    } else if (length(area_lev) == 0) {
      # otherwise, use most common area level in survey_circumcision
      area_lev_df <- survey_circumcision
    }
    
    # pull area_level from area_id col (area_level column not guaranteed)
    area_levs <- area_lev_df %>% 
      dplyr::filter(.data$iso3 == cntry) %>% 
      dplyr::select(contains("area_id")) %>% 
      pull() %>% 
      substr(5, 5) # i.e. LSO_1_01 gives area_level 1
    
    # tabulate and return most common level (i.e. with max count in table)
    area_levs_tbl <- table(area_levs)
    area_lev <- as.numeric(
      names(area_levs_tbl)[area_levs_tbl == max(area_levs_tbl)]
    )
    
    return(area_lev)
  }
  
  # if survey_circumcision is a list or contains more than one country,
  # apply function recursively for each iso3
  is_list <- inherits(survey_circumcision, "list")
  
  # Check if cluster & individuals data is provided; if not, assume sufficient
  # data is included in survey_circumcision alone
  is_add_data_present <- !is.null(survey_clusters) &
    !is.null(survey_individuals)
  
  # split based on iso3 if not already a list and containing > 1 country
  if (!is_list && length(unique(survey_circumcision$iso3)) != 1) {
    survey_circumcision <- split(survey_circumcision, survey_circumcision$iso3)
    
    # arrange and filter for country to ensure splitting is the same
    equal_splits <- function(x, survey_circumcision) {
      x <- x[names(x) %chin% names(survey_circumcision)] # same names
      return(x[order(names(survey_circumcision))]) # order names
    }
    
    # If cluster & individuals data are null, splitting would produce an error
    if (is_add_data_present) {
      survey_individuals <- split(survey_individuals, survey_individuals$iso3)
      survey_individuals <- equal_splits(
        survey_individuals, survey_circumcision
      )
      survey_clusters <- split(survey_clusters, survey_clusters$iso3)
      survey_clusters <- equal_splits(survey_clusters, survey_circumcision)
    }
  }
  
  if (inherits(survey_circumcision, "list")) {
    message(paste0(
      "survey_circumcision supplied is a list and/or contains multiple",
      " countries, applying function recursively for each country present"
    ))
    
    # loop over each country
    surveys <- lapply(seq_along(survey_circumcision), function(i) {
      
      # pull country
      cntry <- unique(survey_circumcision[[i]]$iso3)
      
      message("Preparing surveys for ", cntry, "\n")
      
      if (cntry == "LBR" && cens_age > 29) {
        cens_age <- 29
        message("for LBR, age must be censored to those under 29 (at least)")
      }
      
      # pull latest and first censoring year from survey_id
      survey_years <- as.numeric(substr(unique(
        survey_circumcision[[i]]$survey_id
      ), 4, 7))
      cens_year <- max(survey_years)
      start_year <- max(min(survey_years), start_year)
      
      # filter specific country entries in additional dfs, if they are provided
      if (is_add_data_present) {
        survey_individuals_spec <- dplyr::filter(
          survey_individuals[[i]], .data$iso3 == cntry
        )
        survey_clusters_spec <- dplyr::filter(
          survey_clusters[[i]], .data$iso3 == cntry
        )
      } else {
        survey_individuals_spec <- survey_clusters_spec <- NULL
      }
      survey_circumcision_spec <- dplyr::filter(
        survey_circumcision[[i]], .data$iso3 == cntry
      )
      
      # apply function recursively for each country
      prepare_survey_data(
        areas = dplyr::filter(areas, .data$iso3 == cntry),
        survey_circumcision = dplyr::filter(
          survey_circumcision[[i]], .data$iso3 == cntry
        ),
        survey_individuals = survey_individuals_spec,
        survey_clusters = survey_clusters_spec,
        area_lev = area_lev,
        start_year = start_year,
        cens_year = cens_year,
        cens_age = cens_age,
        rm_missing_type = rm_missing_type,
        norm_kisk_weights = norm_kisk_weights,
        strata.norm = strata.norm,
        strata.kish = strata.kish
      )
    })
    names(surveys) <- names(survey_circumcision)
    return(surveys)
  }
  
  # if we set cens_year == TRUE, calculate cens_year as max survey year
  if (!is.null(cens_year) && cens_year == TRUE) {
    message("cens_year set to TRUE, calculated as last survey year")
    if (!is.null(survey_clusters)) {
      surveys <- survey_clusters$survey_id
    } else {
      surveys <- survey_circumcision$survey_id
    }
    cens_year <- max(as.numeric(
      substr(unique(surveys), 4, 7)
    ))
  }
  
  # remove area_level column if present, to avoid row duplication
  survey_circumcision <- dplyr::select(
    survey_circumcision,
    -dplyr::matches("area_level")
  )
  
  
  # Merging circumcision and individuals survey datasets ---------------------

  # pull original surveys
  orig_surveys <- unique(survey_circumcision$survey_id)
  
  # Join survey data sources together, if they been provided as joined already
  if (is_add_data_present) {
    
    # change colnames to those in line with areas
    if ("geoloc_area_id" %chin% names(survey_clusters)) {
      survey_clusters <- survey_clusters %>%
        dplyr::rename(area_id = .data$geoloc_area_id)
    }
    
    # Merging datasets, if required
    survey_circumcision <- survey_circumcision %>%
      # Merging on individual information to  the circumcision dataset
      dplyr::left_join(
        survey_individuals %>%
          dplyr::select(
            .data$survey_id, .data$cluster_id, .data$individual_id,
            .data$sex, .data$age, .data$indweight
          ),
        by = c("survey_id", "individual_id")
      ) %>%
      # Merging on cluster information to the circumcision dataset
      dplyr::left_join(
        (survey_clusters %>%
           dplyr::mutate(area_id = as.character(.data$area_id)) %>%
           dplyr::select(.data$survey_id, .data$cluster_id, .data$area_id)),
        by = c("survey_id", "cluster_id")
      )
  }
  
  # Add column of NAs for missing age columns
  age_cols <- c("circ_age", "age")
  if (sum(age_cols %chin% names(survey_circumcision)) < 2) {
    missing_age_cols <- age_cols[!age_cols %chin% names(survey_circumcision)]
    survey_circumcision[, missing_age_cols] <- NA
  }
  
  # Remove those with missing circumcison status
  survey_circumcision <- survey_circumcision %>%
    dplyr::filter(
      !is.na(.data$circ_status),
      # need at least one age value for each individual to left censor
      !(is.na(.data$circ_age & is.na(.data$age))),
      !is.na(.data$indweight)
    ) %>%
    # Variables needed for analysis
    dplyr::mutate(
      # Survey year
      year = as.numeric(substr(.data$survey_id, 4, 7)),
      # Year of Birth (estimated as no DOB filly yet)
      yob = .data$year - .data$age,
      # If circumcision age > age of the individual set, reset circumcision age
      circ_age = ifelse(.data$circ_age > .data$age, NA_real_, .data$circ_age)
    )
  
  # Censoring if necessary ---------------------------------------------------
  
  # Censoring at cens_year if assumed no circumcisions after a certain year
  if (!is.null(cens_age)) {
    survey_circumcision <- survey_circumcision %>%
      # Censoring individuals from analysis at cens_age
      dplyr::mutate(
        # No circumcision after cens_age
        circ_status = ifelse(.data$circ_status == 1 &
                               !is.na(.data$circ_age) &
                               .data$circ_age > cens_age, 0, .data$circ_status),
        # Resetting age at circumcision
        circ_age = ifelse(.data$circ_age > cens_age, NA,
                          .data$circ_age
        ),
        # Resetting age for everyone else
        age = ifelse(.data$age > cens_age, cens_age,
                     .data$age
        ),
        # Year of circ/censoring (estimated using the age as no date of circ)
        yoc = ifelse(!is.na(.data$circ_age), .data$yob + .data$circ_age,
                     .data$yob + .data$age
        )
      )
  }
  
  # Censoring at cens_year if assumed no circumcisions after a certain year
  if (!is.null(cens_year)) {
    survey_circumcision <- survey_circumcision %>%
      # Censoring at cens_year
      dplyr::filter(.data$yob < cens_year) %>%
      # Final variables for modelling
      dplyr::mutate(
        # Censoring circumcision status for those circumcised in cens_year,
        # Assuming interval censored people were circumcised before cens_year
        circ_status = ifelse(.data$yoc >= cens_year &
          .data$circ_status == 1 & !is.na(.data$circ_age),
        0.0, .data$circ_status
        ),
        # circ censoring year / censor year in cens_year - 1 at cens_year - 1
        yoc = ifelse(.data$yoc == cens_year, cens_year - 1, .data$yoc)
      )
  }

  # Setting desired level aggregation ----------------------------------------

  if (inherits(areas, "sf")) areas <- sf::st_drop_geometry(areas)
  areas <- dplyr::select(
    areas, .data$area_id, .data$area_name,
    .data$parent_area_id, .data$area_level
  )

  # Getting the area level id to province
  area_lev <- select_area_lev(
    area_lev, 
    cntry = survey_circumcision$iso3[1], 
    is_add_data_present      
  )

  for (i in seq_len(max(areas$area_level))) {
    survey_circumcision <- survey_circumcision %>%
      ## Merging on boundary information
      dplyr::select(-dplyr::matches("area_name")) %>%
      dplyr::left_join(areas, by = "area_id") %>%
      ## Altering area
      dplyr::mutate(
        area_id = dplyr::if_else(
          .data$area_level == area_lev,
          as.character(.data$area_id),
          as.character(.data$parent_area_id)
        )
      ) %>%
      dplyr::select(
        -c(.data$parent_area_id, .data$area_name, .data$area_level)
      )
  }

  # Final preparation of circumcision variables ------------------------------

  # Preparing circumcision variables for the model
  survey_circumcision <- survey_circumcision %>%
    # Merging on the region index
    # Note: inner_join will remove observations that don't
    #       corresponding to a location in "areas"
    dplyr::inner_join(
      dplyr::select(areas, .data$area_id, .data$area_name, .data$area_level),
      by = "area_id"
    ) %>%
    dplyr::mutate(
      # Time interval for the individual
      time1 = .data$yob - start_year + 1,
      time2 = .data$yoc - start_year + 1,
      # Event type
      event = ifelse(
        # event 1 == circumcised and circumcision age known, 0 uncircumcised
        .data$circ_status == 1 & !is.na(.data$circ_age), 1,
        # event 2 == left censored, event 0 == uncircumcised/censored
        ifelse((.data$circ_status == 1 & is.na(.data$circ_age)), 2, 0)
      ),
      # Circumcision age
      circ_age = .data$yoc - .data$yob,
      age = .data$circ_age + 1
    )
  
  # give message on censored individuals
  event_tbl <- table(survey_circumcision$event)
  names(event_tbl) <- c(
    "Uncircumised/right censored (event 0)",
    "Uncensored (event 1)",
    "Left censored (event 2)"
  )[as.numeric(names(event_tbl)) + 1]
  message(
    paste0(utils::capture.output(event_tbl), collapse = "\n")
  )
  
  # Adding circumcision type to dataset
  survey_circumcision <- survey_circumcision %>%
    # Type of circumcision
    dplyr::mutate(
      circ_who = ifelse(.data$circ_who == "other",
                        NA_character_,
                        .data$circ_who
      ),
      circ_where = ifelse(.data$circ_where == "other",
                          NA_character_,
                          .data$circ_where
      ),
      type = dplyr::case_when(
        .data$circ_who == "medical" | .data$circ_where == "medical" ~ "MMC",
        .data$circ_who == "traditional" |
          .data$circ_where == "traditional" ~ "TMC",
        TRUE ~ "Missing"
      )
    )
  
  # Get surveys without any type information
  if (rm_missing_type == TRUE) {
    tmp <- with(
      survey_circumcision, 
      as.data.frame(table(survey_id, type)) # tabulate type
    ) %>%
      # calculate percentage and find surveys with all missing data
      dplyr::group_by(.data$survey_id) %>%
      dplyr::mutate(Freq = .data$Freq / sum(.data$Freq)) %>%
      dplyr::filter(.data$type == "Missing", .data$Freq == 1)
    
    # return message detailing all surveys which are missing
    n <- nrow(tmp)
    if (n > 0) {
      survey_id <- tmp$survey_id
      if (n == 1) {
        message(paste(
          survey_id[1], "has all type == \"Missing\"",
          "and will be removed"
        ))
      } else {
        message(
          paste0(
            paste(paste(survey_id[1:(n - 1)], collapse = ", "),
                  survey_id[n],
                  sep = " & "
            ),
            " have all type == \"Missing\", and will be removed"
          )
        )
      }
    }
    # Removing surveys and individuals without any type information
    survey_circumcision <- survey_circumcision %>%
      dplyr::filter(
        !(.data$survey_id %chin% !!as.character(tmp$survey_id)),
        !(.data$circ_status == 1 & .data$type == "Missing")
      )
  }
  
  # normalise survey weights and apply Kish coefficients, if desired
  if (norm_kisk_weights) {
    survey_circumcision <- normalise_weights_kish(
      survey_circumcision,
      strata.norm = strata.norm,
      strata.kish = strata.kish
    )
  }
  
  # return message for which surveys are discarded & kept (if any)
  remaining_surveys <- unique(survey_circumcision$survey_id)
  if (length(remaining_surveys) != length(orig_surveys)) {
    removed_surveys <- orig_surveys[!orig_surveys %chin% remaining_surveys]
    surveys <- list(removed_surveys, remaining_surveys)
    lengths <- c(length(removed_surveys), length(remaining_surveys))
    initial <- c("Surveys removed: ", "Surveys remaining: ")
    invisible(lapply(seq_along(surveys), function(i) {
      if (lengths[i] == 1) {
        message(paste0(initial[i], surveys[[i]]))
      } else {
        message(
          paste0(
            initial[i],
            paste(paste(
              surveys[[i]][1:(lengths[i] - 1)],
              collapse = ", "
            ),
            surveys[[i]][lengths[i]],
            sep = " & "
            )
          )
        )
      }
    }))
  }
  
  # Returning prepped circumcision datasets
  return(survey_circumcision)
}
