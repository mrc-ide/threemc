#' @title PSNU Area Levels for SSA
#' @description PSNU area levels for Sub-Saharan African countries. These
#' are the recommended levels at which to perform modelling etc., for each
#' respective country. Inferences on larger regions (i.e. lower PSNU area
#' levels) can be made by aggregating results for higher area levels. The
#' dataset contains the following fields:
#' \itemize{
#'   \item{\code{iso3}}{character ISO3 codes for Sub-Saharan African
#'   countries.}
#'   \item{\code{psnu_area_level}}{integer The sub national level considered
#'   to be the organizational level in which a country has prioritised their
#'   program. Increasing values refer to more granular regional distinctions.}
#' }
#' @format A \code{data.frame} with 29 rows and 2 variables:
#' @docType data
#' @usage data(datapack_psnu_area_level)
#' @keywords datasets
#' @name datapack_psnu_area_level
NULL

#' @title WCA - ESA key for Sub-Saharan African countries
#' @description Western and Central Africa (WCA) - Eastern and Southern Africa 
#' (ESA) categorisation for Sub-Saharan African countries. Also includes 
#' North-South-East-West categorisation. 
#' \itemize{
#'   \item{\code{iso3}}{character ISO3 codes for Sub-Saharan African
#'   countries.}
#'   \item{\code{region}}{character ESA-WCA categorisation for each 
#'   \code{iso3}}
#'   \item{\code{four_region}}{character North-South-East-West categorisation 
#'   for each \code{iso3}}
#' }
#' @format A \code{data.frame} with 38 rows and 3 variables:
#' @docType data
#' @usage data(esa_wca_regions)
#' @keywords datasets
#' @name esa_wca_regions
NULL
