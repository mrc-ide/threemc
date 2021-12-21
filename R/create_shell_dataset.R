#' @title Create Shell Dataset for Estimating Empirical Circumcision Rate
#'
#' @description Create a matrix to estimate the cumulative hazard rate needed
#' for survival analysis by age and time. The option to include an additional 
#' stratification variable is also available, creating a 3D hazard function.
#' 
#' @param dat Dataset used for modelling.
#' @param subset Subset for dataset.
#' @param time1 Variable name for time of birth.
#' @param time2 Variable name for time circumcised or censored.
#' @param timecaps Window to fix temporal dimension before and after.
#' @param Ntime Number of time points (if NULL, function will calculate).
#' @param age - Variable with age circumcisied or censored.
#' @param Nage Number of age groups (if NULL, function will calculate).
#' @param strat Variable to stratify by in using a 3D hazard function.
#' @param Nstrat Number of stratification groups (if NULL, function will 
#' calculate).
#' 
#' @return Matrix for selecting instananeous hazard rate.
#' @export
#' 
#' @import dplyr
#' @importFrom tidyr crossing

# function to create shell dataset for estimating the empirical circ rate
create_shell_dataset <- function(survey_circumcision, 
                                 areas, 
                                 area_lev,
                                 time1  = "time1",
                                 time2  = "time2",
                                 strat  = "space",
                                 age    = "age",
                                 circ   = "indweight_st") {
  
  #' !! JE: the area_id field here used the area_id that appeared in the
  #'        circumcision dataset. If there were some area_id with no
  #'        observations, then they were dropped, which created a
  #'        misalignment of the indexing.  Instead, use the areas dataset
  #'        to construct this output frame to ensure all districts are
  #'        represented.
  #'
  #'        I suspect that we also want the circ_age field to be constructed
  #'        based on the theoretical maximum circumcision age that we want
  #'        outputs for, rather than the maximum observed age; but not 100%
  #'        sure.
  
  
  if (missing(area_lev)) {
    warning("area_lev arg missing, taken as maximum area level in areas")
    area_lev <- max(areas$area_level, na.rm = T)
  }
  # remove spatial elements from areas, take only specified/highest area level
  if ("sf" %in% class(areas)) {
    areas_model <- st_drop_geometry(areas)
  } else areas_model <- areas
  areas_model <- areas_model %>%
    filter(area_level == area_lev) %>%
    select(area_id, area_name, space)
  
  #' create skeleton dataset with row for every unique area_id, area_name, 
  #' space, year and circ_age
  out <- crossing(areas_model,
                  year = seq(2006, 2021, by =  1),
                  circ_age = 0 : max(survey_circumcision$circ_age)) %>%
    # Getting time and age variable
    mutate(time = year - 2006 + 1,
           age  = circ_age + 1) %>%
    # Sorting dataset
    arrange(space, age, time)
  
  # Obtain N person years
  out_int_mat <- create_integration_matrix_agetime(
    dat = survey_circumcision,
    time1 = time1,
    time2 = time2,
    strat = strat,
    age   = age,
    Ntime = length(unique(out$time))
  )
  out$N <- as.vector(survey_circumcision$indweight_st %*% out_int_mat)
  
  #' calculate empirical agetime hazard matrices for different circumcision 
  #' types, and take column sums:
  subsets <- c("event == 1 & type == 'MMC'", # N MMC 
               "event == 1 & type == 'TMC'", # N TMC
               "event == 0", # N censored (i.e. not circumcised)
               "event == 2" # N left-censored (circ at unknown age)
  )
  agetime_hazard_matrices <- lapply(subsets, function(x) {
    create_hazard_matrix_agetime(
      dat = survey_circumcision,
      subset = x,
      time1 = time1,
      time2 = time2,
      strat = strat,
      age   = age,
      circ  = circ,
      Ntime = length(unique(out$time))
    ) %>%
      colSums()
  })
  
  # add to out:
  empirical_circ_cols <- c("obs_mmc", "obs_tmc", "cens", "icens")
  out[, empirical_circ_cols] <- agetime_hazard_matrices
  
  return(out)
}
