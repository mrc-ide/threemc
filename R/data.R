#' @title PSNU Area Levels for SSA
#' @description PSNU area levels for Sub-Saharan African countries. These
#' are the recommended levels at which to perform modelling etc., for each
#' respective country. Inferences on larger regions (i.e. lower PSNU area
#' levels) can be made by aggregating results for higher area levels. The
#' dataset contains the following fields:
#' \itemize{
#'   \item{\code{iso3}}{character ISO3 codes for various Sub-Saharan African
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
