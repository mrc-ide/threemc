
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
#' @param probs Specific quantiles to be calculated,
#' Default: c(0.5, 0.025, 0.975)
#' @param names Parameter with \link[stats]{quantile}: logical; if true, the
#' result has a names attribute. Set to FALSE for speedup with many probs,
#' Default: FALSE
#' @param ... Further arguments passed to or from other methods.
#' @return Input `out` `data.frame`, including columns with quantiles for
#' hazard rates etc for different circumcision types, and for overall
#' circumcision.

#' @seealso
#'  \code{\link[threemc]{circ_sample_tmb}}
#'  \code{\link[stats]{quantile}}
#' @rdname compute_quantiles
#' @export

## function to calculate quantiles for rates and cumulative hazard
## maybe add "textProgressBar" to for loop if it's too slow?
compute_quantiles <- function(out,
                              fit,
                              probs = c(0.5, 0.025, 0.975),
                              names = FALSE, ...) {
  probs <- sort(probs)

  ## function to do "legwork" of computing quantiles
  quantile_fun <- function(.data, probs = probs, names = names, ...) {
    .data <- t(apply(.data, 1, quantile, probs = probs, names = names, ...))
  }

  ## want to add quantile columns to out for following hazards
  types <- c(
    "rate_mmc", "rate_tmc", "rate", "surv", # rates
    "inc_tmc", "inc_mmc", "inc", # incidence
    "cum_inc_tmc", "cum_inc_mmc", "cum_inc"
  ) # cumulative incidence
  ## append "M" (median), "L" (lower) and "U" (upper) to these column names
  types <- lapply(types, function(x) paste0(rep(x, 3), c("M", "L", "U")))

  ## pull samples for each of these hazards from fit object
  samples <- with(
    fit$sample,
    list(
      haz_mmc, haz_tmc, haz_mmc + haz_tmc, surv,
      inc_tmc, inc_mmc, inc_tmc + inc_mmc,
      cum_inc_tmc, cum_inc_mmc, cum_inc_tmc + cum_inc_mmc
    )
  )
  ## for each hazard, calculate quantiles and add as appropriate cols to out
  for (i in seq_along(types)) {
    out[, c(types[[i]])] <- quantile_fun(samples[[i]],
      probs = probs,
      names = names, ...
    )
  }

  return(out)
}
