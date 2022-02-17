#' @title Sample TMB fit for Circumcision Model
#' @description  Sample from TMB object, using \link[naomi]{sample_tmb}. Saves
#' changing object to "Naomi" format. Also produces and returns standard
#' deviation report outputted by \link[TMB]{sdreport}.
#'
#' @param obj TMB object/AD model outputted by \link[TMB]{MakeADFun}.
#' @param opt Optimised TMB model, outputted by optimisation function such
#' as \link[stats]{nlminb} or \link[stats]{optim}.
#' @param nsample Number of samples to be generated, Default: 1000
#' @param ...  Further arguments passed to \link[naomi]{sample_tmb}.
#' @return Object of class "naomi_fit", containing the original TMB object
#' ("obj"), the standard deviation report for optimised AD model (from
#' \link[TMB]{sdreport}) and `n_samples` samples for the (cumulative) incidence
#' and hazard rate of circumcision for the region(s) in question.
#'
#' @seealso
#'  \code{\link[TMB]{sdreport}}
#'  \code{\link[naomi]{sample_tmb}}
#' @rdname circ_sample_tmb
#' @export
circ_sample_tmb <- function(obj, opt, nsample = 1000, ...) {

  ## Getting the TMB into "Naomi" format to sample from using the NAOMI package
  opt$par.fixed <- opt$par
  opt$par.full <- obj$env$last.par
  fit <- c(opt, obj = list(obj))
  class(fit) <- "naomi_fit"

  ## Look at standard deviation report
  fit$sdreport <- TMB::sdreport(fit$obj, fit$par, getJointPrecision = TRUE)

  ## Generating samples
  fit <- naomi::sample_tmb(fit, nsample = nsample, ...)
  
  # ensure names for MC columns in fit have the suffix "_mc"
  fit$sample <- append_mc_name(fit$sample)
  
  return(fit)
}
