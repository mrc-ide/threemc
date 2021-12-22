#' @title Sample TMB fit for Circumcision Model
#' @description  Sample from TMB object, using \link[naomi]{sample_tmb}. Saves
#' changing object to "Naomi" format. Also produces and returns standard
#' deviation report outputted by \link[TMB]{sdreport}.
#'
#' @param obj TMB object/AD model outputted by \link[TMB]{MakeADFun}.
#' @param opt Optimised TMB model, outputted by optimisation function such
#' as \link[stats]{nlminb} or \link[stats]{optim}.
#' @param n_samples Number of samples to be generated. Additional parameters
#' for \link[naomi]{sample_tmb} can also be supplied, Default: 1000
#' @return Object of class "naomi_fit", containing the original TMB object
#' ("obj"), the standard deviation report for optimised AD model ("sdreport")
#' and `n_samples` samples for the (cumulative) incidence and hazard rate of
#' circumcision for the region(s) in question.
#'
#' @seealso
#'  \code{\link[TMB]{sdreport}}
#'  \code{\link[naomi]{sample_tmb}}
#' @rdname circ_sample_tmb
#' @export
#'
#' @importFrom TMB sdreport
#' @importFrom naomi sample_tmb
circ_sample_tmb <- function(obj, opt, n_samples = 1000, ...) {
  
  # Getting the TMB into "Naomi" format to sample from using the NAOMI package
  opt$par.fixed <- opt$par
  opt$par.full <- obj$env$last.par
  fit <- c(opt, obj = list(obj))
  class(fit) <- "naomi_fit"

  # Look at standard deviation report
  fit$sdreport <- sdreport(fit$obj, fit$par, getJointPrecision = TRUE)
  
  # Generating samples
  fit <- sample_tmb(fit, n_samples = 1000, ...)
  return(fit)
}
