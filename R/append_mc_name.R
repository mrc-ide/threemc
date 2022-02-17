#' @title Append "mc" to Appropriate Column names
#' @description Ensure names for MC columns in fit have the suffix "_mc"
#' @param .data Dataframe/tibble whose columns include "MC" calculations.
#' @return \code{.data}, with column names appended appropriately.
#' @export
append_mc_name <- function(.data) {
  mmc_tmc <- paste(c("mmc", "tmc"), collapse = "|")
  locs <- !(grepl(paste(c(mmc_tmc, "mc"), collapse = "|"), names(.data)))
  names(.data)[locs] <- paste0(names(.data)[locs], "_mc")

  return(.data)
}
