#' @title PSNU Area Levels for SSA 
#' @description PSNU area levels for Sub-Saharan African countries. These
#' are the recommended levels at which to perform modelling etc., for each 
#' respective country. Inferences on larger regions (i.e. lower PSNU area 
#' levels) can be made by aggregating results for higher area levels.
#' @format A data frame with 29 rows and 2 variables:
#' \describe{
#'   \item{\code{iso3}}{character ISO3 codes for various Sub-Saharan African 
#'   countries.}
#'   \item{\code{psnu_area_level}}{integer The sub national level considered 
#' to be the organizational level in which a country has prioritised their 
#' program. Increasing values refer to more granular regional distinctions.} 
#'}
"datapack_psnu_area_level"