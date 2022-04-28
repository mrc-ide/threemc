#' @title Produce Population weighted Aggregated Samples
#' @description Aggregate by area, year, age and type (weighted by population),
#' and convert to a percentage/probability.
#' @param .data \code{data.frame} of unaggregated modelling results. 
#' @param fit \code{TMB} list containing model parameters, nested list of 
#' `samples` for the (cumulative) incidence and hazard rate of circumcision for 
#' the region(s) in question.
#' @param areas `sf` shapefiles for specific country/region. 
#' @param populations \code{data.frame} containing populations for each
#' region in tmb fits.
#' @param age_var Determines whether you wish to aggregate by discrete ages or 
#' age groups (0-4, 5-9, 10-14, and so on). 
#' @param type Determines which aspect of MC in the regions in question you wish
#' to aggregate for. Can be one of "probability", "incidence" or "prevalence".
#' @param area_lev PSNU area level for specific country. 
#' @param N Number of samples to be generated, Default: 100
#' @param prev_year If type == "prevalence", choose year to compare prevalence 
#' with. 
#' @return \code{data.frame} with samples aggregated by \code{aggr_cols} and
#' weighted by population.
#' @importFrom dplyr %>% 
#' @importFrom rlang .data
#' @rdname threemc_aggregate
#' @export
threemc_aggregate <- function(
  # datasets
  .data, fit, areas, populations, # datasets
  # options
  age_var = c("age", "age_group"), 
  type = c("probability", "incidence", "prevalence"), area_lev, N = 100, 
  prev_year = 2008
  ) {
  
  #### Preparing location/shapefile information ####
  if (inherits(areas, "sf")) {
    areas <- sf::st_drop_geometry(areas)
  }
  
  areas <- areas %>%
    # Add a unique identifier within Admin code and merging to boundaries
    dplyr::group_by(.data$area_level) %>%
    dplyr::mutate(space = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # wide formatted areas, for changing area levels later
  areas_wide <- areas %>%
    dplyr::select(area_id, area_name, parent_area_id, area_level) %>%
    naomi::spread_areas()
  
  # Model with Probability of MC
  .data$model <- "No program data"
  fit_no_prog <- fit # need to load model with programme data as well
  
  #### Load rates from survival model ####
  .data <- prepare_sample_data(N = N,
                                 populations = populations,
                                 no_prog_.data = .data,
                                 no_prog_tmb_fit = fit_no_prog,
                                 type = type)
  
  rm(fit, fit_no_prog, populations); gc()
  
  #### Prepare .data for output ####
  
  # collect .data for lower area hierarchies by joining higher area
  # hierarchies (do so until you reach "area 0")
  .data <- combine_areas(.data, areas_wide, area_lev, join = FALSE)
  
  # aggregate samples for each individual age or age group
  if (age_var == "age") {
    .data <- aggregate_sample(.data)
  } else .data <- aggregate_sample_age_group(.data)
  
  # additional aggregations to perform for prevalence
  if (type == "prevalence") {
    # calculate change in prevalence since prev_year
    data_change_prev_year <- prevalence_change(
      .data, spec_year = prev_year
    )
    
    # Getting number of people circumcised
    data_n <- n_circumcised(.data)
    
    # add "change from prev_year" and "n_circumcised" .data to .data
    .data <- rbind(.data, data_change_prev_year, data_n)
  }
  # calculate summary statistics (mean, sd, quantiles)
  .data <- posterior_summary_fun(.data)
  
  # Merge regional information on the dataset (i.e. parent area info) & return
  .data <- merge_area_info(.data, areas)
}
