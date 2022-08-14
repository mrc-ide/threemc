
#### Main Function #### 

#' @title Produce TMB model fit with sample, or re-sample from existing 
#' optimised model fit. 
#' @description Optimises threemc objective function and produces samples from
#' model fit (if so desired). If provided with an existing optimised model 
#' `fit`, can also perform re-sampling. 
#' @param fit Optional "small" fit object with no `sample`. Specifying `fit` 
#' means you do not need to specify `dat_tmb` or `parameters`, as argument 
#' specifications will be overridden by those stored in `fit`. 
#' @param dat_tmb \code{list} of data required for model fitting, outputted
#' by \link[threemc]{threemc_prepare_model_data}, which includes:
#' \itemize{
#'   \item{design_matrices}{Includes`X_fixed_mmc`, `X_fixed_tmc`, `X_time_mmc`,
#'   `X_age_mmc`, `X_age_tmc`, `X_space_mmc`, `X_space_tmc`, `X_agetime_mmc`,
#'   `X_agespace_mmc`, `X_agespace_tmc`, `X_spacetime_mmc`. Design
#'    Create design matrices for fixed effects and temporal, age, space and
#'    interaction random effects}
#'    \item{integration matrices}{Includes `IntMat1`, `IntMat2`. Integration
#'    matrices for selecting the instantaneous hazard rate.}
#'    \item{survival matrices}{Includes `A_mmc`, `A_tmc`, `A_mc`, `B`, `C`.
#'    Survival matrices for MMC, TMC, censored and left censored}
#'    \item{Q_space}{Precision/Adjacency matrix for the spatial random effects.
#'    }
#' }
#' @param parameters \code{list} of fixed and random model parameters.
#' @param maps \code{list} of factors with value NA, the names of which
#' indicate parameters to be kept fixed at their initial value throughout the
#' optimisation process.
#' @param randoms \code{vector} of random effects.
#' @param mod TMB model, one of either
#' "Surv_SpaceAgeTime_ByType_withUnknownType" or "Surv_SpaceAgeTime" if the
#' surveys for the country in question make no distinction between circumcision
#' type (i.e whether they were performed in a medical or traditional setting).
#' @param sample If set to TRUE, has function also return N samples for
#' medical, traditional and total circumcisions, Default: TRUE
#' @param smaller_fit_obj Returns a smaller fit object. Useful for saving the 
#' fit object for later aggregations. 
#' @param sdreport If set to TRUE, produces the standard deviation report for
#' the model, Default: FALSE
#' @param N Number of samples to be generated, Default: 1000
#' @param ... Further arguments passed to internal functions.
#' @return TMB model fit, including optimised parameters, hessian matrix, 
#' samples and standard deviation report (if desired).
#' @rdname threemc_fit_model
#' @export
threemc_fit_model <- function(
  fit = NULL, dat_tmb = NULL, mod, parameters = NULL, maps = NULL,
  randoms = c(
    "u_time_mmc", "u_age_mmc", "u_space_mmc",
    "u_agetime_mmc", "u_agespace_mmc", "u_spacetime_mmc",
    "u_age_tmc", "u_space_tmc", "u_agespace_tmc"
  ),
  sample = TRUE, smaller_fit_obj = FALSE, sdreport = FALSE, N = 1000, ...
) {
  
  # for specified "smaller fit" object (i.e. fit which requires resampling)
  if (!is.null(fit)) {
    if (!is.null(fit$sample)) stop("Sample already present in fit object")
    if (!is.null(dat_tmb) | !is.null(parameters)) {
      message(paste0(
       "No need to specify dat_tmb or parameters for non-null fit, as they are",
       " replaced by those stored in fit"
      ))
   }
   # pull dat_tmb and parameters from small fit 
   dat_tmb <- fit$tmb_data
   parameters <- split(fit$par.full, names(fit$par.full))
   init_params <- fit$par_init
   # pull different parameters depending on whether the model has mmc/tmc split
   if (mod == "Surv_SpaceAgeTime_ByType_withUnknownType") {
     parameters <- parameters[names(fit$par_init)]
   } else {
     # only need names and lengths, not values
     names(init_params) <- stringr::str_remove_all(
       names(init_params), "_mmc|_tmc"
     )
     init_params <- init_params[!duplicated(names(init_params))]
     parameters <- parameters[names(init_params)]
     # remove duplicate parameters
     parameters <- parameters[!duplicated(names(parameters))]
   }
   
   if (!is.null(maps)) {
     # ensure mapped parameters are in the same order as parameters for model
     mapped_pars <- is.na(names(parameters))
     param_order <- names(init_params)[mapped_pars]
     maps <- maps[match(names(maps), param_order)]
     
     # replace NAs in parameters with mapped parameters in par_init
     parameters[mapped_pars] <- init_params[
       names(init_params) %in% names(maps)
     ]
     names(parameters)[mapped_pars] <- names(maps)
   }
   
   is_matrix <- sapply(init_params, is.matrix)
   parameters[is_matrix] <- Map(matrix,
                             parameters[is_matrix],
                             nrow = lapply(init_params[is_matrix], nrow),
                             ncol = lapply(init_params[is_matrix], ncol))
  # if no fit == NULL, must have non-null dat_tmb & parameters
  } else {
    
    if (is.null(dat_tmb) | is.null(parameters)) {
      
      stop("Please specify non-null dat_tmb and parameters")
    }
  }

  # remove "mmc" from parameter & matrix names if required
  if (mod == "Surv_SpaceAgeTime") {
    remove_type_distinction <- function(x) {
      names(x) <- stringr::str_remove(names(x), "_mmc")
      x <- x[!grepl("_tmc", names(x))]
    }

    dat_tmb <- remove_type_distinction(
      dat_tmb[!names(dat_tmb) %in% c("A_mmc", "A_tmc")]
    )
    names(dat_tmb)[names(dat_tmb) == "A_mc"] <- "A"

    parameters <- remove_type_distinction(parameters)
    randoms <- unique(stringr::str_remove(randoms, "_tmc|_mmc"))
  }

  # Only have named random parameters
  randoms <- randoms[randoms %in% names(parameters)]
  if (length(randoms) == 0) randoms <- NULL

  # Create TMB object
  obj <- TMB::MakeADFun(
    dat_tmb,
    parameters,
    random = randoms,
    map = maps,
    method = "BFGS",
    hessian = TRUE,
    DLL = mod,
    ...
  )
  # for specified fit, simply resample and return
  if (!is.null(fit)) {
    fit$obj <- obj
    fit$obj$fn()
    fit <- circ_sample_tmb(
      fit = fit, obj = obj, nsample = N, sdreport = sdreport
    )
    fit$tmb_data <- fit$par_init <- NULL # make fit object smaller for saving
    return(fit)
  }

  # Run optimiser (use optim if all pars are fixed, nlminb otherwise)
  if (length(obj$par) == 0) {
    opt <- do.call(stats::optim, obj, ...)
  } else {
    opt <- stats::nlminb(
      start   = obj$par,
      obj     = obj$fn,
      gr      = obj$gr,
      control = list(trace = 1),
      ...
    )
  }
  
  # sample from TMB fit
  if (sample == TRUE) {
    fit <- circ_sample_tmb(
      obj = obj, opt = opt, nsample = N, sdreport = sdreport
    )
    # return smaller fit object
    if (smaller_fit_obj == TRUE) {
      fit <- minimise_fit_obj(fit, dat_tmb, parameters)
    }
    return(fit)
  } else {
    return(opt)
  }
}

#### circ_sample_tmb ####

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
#' @keywords internal
circ_sample_tmb <- function(
  fit = NULL, obj = NULL, opt, sdreport = FALSE, nsample = 1000, ...
  ) {

  # Getting the TMB into "Naomi" format to sample from using the NAOMI package
  if (is.null(fit)) {
    opt$par.fixed <- opt$par
    opt$par.full <- obj$env$last.par
    fit <- c(opt, obj = list(obj))
  } 
  class(fit) <- "naomi_fit"  
  
  # Look at standard deviation report
  if (sdreport == TRUE) {
      fit$sdreport <- TMB::sdreport(fit$obj, fit$par, getJointPrecision = TRUE)
  }

  # Generating samples
  fit <- naomi::sample_tmb(fit, nsample = nsample, ...)

  # ensure names for MC columns in fit have the suffix "_mc"
  fit$sample <- append_mc_name(fit$sample)

  return(fit)
}

#### minimise_fit_obj #### 

#' @title Minimise Fit Object Size
#' @description Return minimised fit object. Often useful when saving the fit 
#' object for later aggregation. 
#' @param fit Fit object returned by \link[naomi]{sample_tmb}, which includes,
#' among other things, the optimised parameters and subsequent sample for our 
#' TMB model.
##' @param dat_tmb \code{list} of data required for model fitting, outputted
#' by \link[threemc]{threemc_prepare_model_data}.
#' @param parameters \code{list} of fixed and random model parameters.
#' @return Object of class "naomi_fit".
#' @rdname minimise_fit_obj
#' @export
minimise_fit_obj <- function(fit, dat_tmb, parameters) {
  
  fit_small <- fit
  fit_small$tmb_data <- dat_tmb
  fit_small$par_init <- parameters
  fit_small$sample <- NULL
  fit_small$obj <- NULL
  
  return(fit_small)
}
