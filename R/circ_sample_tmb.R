# function to sample from TMB fit
circ_sample_tmb <- function(obj, opt) {
  
  # Getting the TMB into "Naomi" format to sample from using the NAOMI package
  opt$par.fixed <- opt$par
  opt$par.full <- obj$env$last.par
  fit <- c(opt, obj = list(obj))
  class(fit) <- "naomi_fit"
  
  # Look at standard deviation report
  fit$sdreport <- sdreport(fit$obj, fit$par, getJointPrecision = TRUE)
  
  # Generating samples
  fit <- sample_tmb(fit)
  return(fit)
}
