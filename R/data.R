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
#' @format A \code{data.frame} with 38 rows and 3 variables.
#' @docType data
#' @usage data(esa_wca_regions)
#' @keywords datasets
#' @name esa_wca_regions
NULL

#' @title Malawi shapefiles
#' @description `sf` shapefile representation of Malawi, as a multipolygon.
#' \itemize{
#'   \item{\code{iso3}}{character ISO3 codes for Sub-Saharan African
#'   countries.}
#'   \item{\code{area_id}}{Unique ID for each region in MWI. Formatted as 
#'   "County_area_level_ID" (e.g. MWI_3_05 for Mzimba)}
#'   \item{\code{area_name}}{Name of region in question}
#'   \item{\code{parent_area_id} {Unique ID for region's parent region.}}
#'   \item{\code{area_level} {Numeric value denoting area level of 
#'   area, in decreasing granularity.}}
#'   \item{\code{area_level_label} {Translates numeric area level to meaning 
#'   in country in question. For example, in Malawi a region of area level 3 
#'   is a "District".}}
#'   \item{\code{area_sort_order} {Order to sort areas in when plotting, 
#'   roughly equivalent to a `geofacet` grid.}}
#'   \item{\code{center_x} {X coordinate for centre of region's multipolygon}}
#'   \item{\code{center_y} {Y coordinate for centre of region's multipolygon}}
#'   \item{\code{geometry}} {`sfc_MULTIPOLYGON` representation of region's  
#'   spatial geometry}
#' }
#' @format A \code{sf} collecton of 6 features and 9 fields, including a  
#' \code{data.frame} with 387 rows, 10 variables, and a sf 
#' @docType data
#' @usage data(demo_areas)
#' @keywords datasets
#' @name demo_areas
NULL

#' @title Malawi populations 
#' @description Single age, aggregated male populations for each area in 
#' Malawi.
#' \itemize{
#'   \item{\code{iso3}}{character ISO3 codes for Sub-Saharan African
#'   countries.}
#'   \item{\code{area_id}}{Unique ID for each region in MWI. Formatted as 
#'   "County_area_level_ID" (e.g. MWI_3_05 for Mzimba)}
#'   \item{\code{area_level} {Numeric value denoting area level of 
#'   area, in decreasing granularity.}}
#'   \item{\code{area_name}}{Name of region in question}
#'   \item{\code{year}}{Year for population in question}
#'   \item{\code{age}}{Age for population in question}
#'   \item{\code{population} {(Male) Population for each unique 
#'   area-year-age combination.}}
#' }
#' @format A \code{data.frame} with 58806 rows and 7 variables.
#' @docType data
#' @usage data(demo_populations)
#' @keywords datasets
#' @name demo_populations
NULL

#' @title Malawi surveys 
#' @description Circumcision surveys for Malawi. 
#' \itemize{
#'   \item{\code{iso3}}{character ISO3 codes for Sub-Saharan African
#'   countries.}
#'   \item{\code{survey_id}}{Survey id for each record.}
#'   \item{\code{area_id}}{Unique ID for each region in MWI. Formatted as 
#'   "County_area_level_ID" (e.g. MWI_3_05 for Mzimba)}
#'   \item{\code{area_level} {Numeric value denoting area level of 
#'   area, in decreasing granularity.}}
#'   \item{\code{age}{Age at interview.}}
#'   \item{\code{dob_cmc}{CMC (Century Month Code) date of birth of 
#'   individual.)}}
#'   \item{\code{interview_cmc}{CMC date of interview.}}
#'   \item{\code{indweight}{Weighting for survey record in question}}
#'   \item{\code{circ_status}{Circumcision status of individual, 1 indicating 
#'   circumcision and 0 indicating right-censoring.}} 
#'   \item{\code{circ_age}{Age at circumcision, if applicable.}}
#'   \item{\code{circ_who}{Circumcision provider, either medical or 
#'   traditional.}}
#'   \item{\code{circ_where}{Circumcision location, either medical or 
#'   traditional.}}
#' }
#' @format A \code{data.frame} with 29313 rows and 12 variables.
#' @docType data
#' @usage data(demo_survey_circumcision)
#' @keywords datasets
#' @name demo_survey_circumcision
NULL
