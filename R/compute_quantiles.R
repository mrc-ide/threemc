
#' @title Calculate Quantiles for Rates and Cumulative Hazard
#' @description Calculate quantiles for samples of rates and cumulative hazard
#' outputted from \link[threemc]{circ_sample_tmb}, and add them as columns to
#' the shell `data.frame` `out` with estimated empirical circumcision rates.
#'
#' @param out Shell dataset with a row for every unique record in
#' circumcision survey data for a given area. Also includes empirical estimates
#'  for circumcision estimates for each unique record.
#' @param fit Object containing `samples` list for the (cumulative) incidence
#' and hazard rate of circumcision for the region(s) in question.
#' @param area_lev PSNU area level for specific country
#' @param probs Specific quantiles to be calculated,
#' Default: c(0.025, 0.5, 0.975)
#' @param names Parameter with \link[stats]{quantile}: logical; if true, the
#' result has a names attribute. Set to FALSE for speedup with many probs,
#' Default: FALSE
#' @param ... Further arguments passed to \link[stats]{quantile}.
#' @return Input `out` `data.frame`, including columns with quantiles for
#' hazard rates etc for different circumcision types, and for overall
#' circumcision.

#' @seealso
#'  \code{\link[stats]{quantile}}
#'  @importFrom dplyr %>%
#'  @importFrom rlang .data
#' @rdname compute_quantiles
#' @export

## function to calculate quantiles for rates and cumulative hazard
## maybe add "textProgressBar" to for loop if it's too slow?
compute_quantiles <- function(out,
                              fit,
                              area_lev = NULL,
                              probs = c(0.025, 0.5, 0.975),
                              names = FALSE,
                              ...) {
  if (length(probs) != 3) stop("probs should be of length 3")
  # ensure that probs are sorted, regardless of input order
  probs <- sort(probs)
  
  # subset to area level of interest, if desired
  if (!is.null(area_lev)) {
    out <- out %>% 
      dplyr::filter(.data$area_level == area_lev)
  }

  # function to do "legwork" of computing quantiles
  quantile_fun <- function(.data, probs = probs, names = names, ...) {
    .data <- t(apply(.data, 1, stats::quantile,
      probs = probs, names = names, ...
    ))
  }

  # want to add quantile columns to out for following hazards
  types <- c(
    "rate_mmc", "rate_tmc", "rate", "surv", # rates
    "inc_mmc", "inc_tmc", "inc", # incidence
    "cum_inc_mmc", "cum_inc_tmc", "cum_inc" # cumulative incidence
  )

  # ensure names for MC columns in fit have the suffix "_mc"
  fit$sample <- append_mc_name(fit$sample)

  # if we are modelling only MC coverage, only want non-type specific "types"
  samples <- NULL
  if (!"haz_mmc" %in% names(fit$sample)) {
    types <- types[!grepl(paste(c("mmc", "tmc"), collapse = "|"), types)]
    # pull corresponding samples for each of these hazards from fit object
    samples <- with(fit$sample, list(haz_mc, surv_mc, inc_mc, cum_inc_mc))
  }

  # append "L" (lower), "M" (mean), and "U" (upper) to these column names
  # probs must be sorted to ensure this is in the right order
  types <- lapply(types, function(x) paste0(rep(x, 3), c("L", "M", "U")))

  # pull samples for each of these hazards from fit object
  if (is.null(samples)) { # only do if explicitly modelling MMC and TMC
    samples <- with(
      fit$sample,
      list(
        haz_mmc, haz_tmc, haz_mc,
        surv_mc,
        inc_mmc, inc_tmc, inc_mc,
        cum_inc_mmc, cum_inc_tmc, cum_inc_mc
      )
    )
  }

  # for each hazard, calculate quantiles and add as appropriate cols to out
  for (i in seq_along(types)) {
    out[, c(types[[i]])] <- quantile_fun(samples[[i]],
      probs = probs,
      names = names,
      ...
    )
  }

  return(out)
}
