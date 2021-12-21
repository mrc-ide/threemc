#' function to calculate quantiles for rates and cumulative hazard
#' maybe add "textProgressBar" to for loop if it's too slow?
compute_quantiles <- function(out, 
                              fit, 
                              probs = c(0.5, 0.025, 0.975), 
                              names = FALSE, ...) {
  
  probs <- sort(probs)
  
  # function to do "legwork" of computing quantiles
  quantile_fun <- function(.data, probs = probs, names = names, ...) {
    .data <- t(apply(.data, 1, quantile, probs = probs, names = names, ...))
  }
  
  # want to add quantile columns to out for following hazards
  types <- c("rate_mmc", "rate_tmc", "rate", "surv", # rates 
             "inc_tmc", "inc_mmc", "inc", # incidence
             "cum_inc_tmc", "cum_inc_mmc", "cum_inc") # cumulative incidence
  # append "M" (median), "L" (lower) and "U" (upper) to these column names
  types <- lapply(types, function(x) paste0(rep(x, 3), c("M", "L", "U"))) 
  
  # pull samples for each of these hazards from fit object
  samples <- with(fit$sample, 
                  list(haz_mmc, haz_tmc, haz_mmc + haz_tmc, surv, 
                       inc_tmc, inc_mmc, inc_tmc + inc_mmc, 
                       cum_inc_tmc, cum_inc_mmc, cum_inc_tmc + cum_inc_mmc)
  )
  # for each hazard, calculate quantiles and add as appropriate cols to out
  for (i in seq_along(types)) {
    out[, c(types[[i]])] <- quantile_fun(samples[[i]], probs = probs, 
                                         names = names, ...)
  }
  
  return(out)
}
