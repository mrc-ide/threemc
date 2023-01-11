
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
#' @param verbose Boolean specifying whether you want detailed updates on 
#' function operations and progress, default = TRUE
#' @param ... Further arguments passed to internal functions.
#' @return TMB model fit, including optimised parameters, hessian matrix,
#' samples and standard deviation report (if desired).
#' @rdname threemc_fit_model
#' @export
threemc_fit_model <- function(fit = NULL,
                              dat_tmb = NULL,
                              mod = NULL,
                              parameters = NULL,
                              maps = NULL,
                              randoms = c(
                                "u_time_mmc", "u_age_mmc", "u_space_mmc",
                                "u_agetime_mmc", "u_agespace_mmc",
                                "u_spacetime_mmc", "u_age_tmc",
                                "u_space_tmc", "u_agespace_tmc"
                              ),
                              sample = TRUE,
                              smaller_fit_obj = FALSE,
                              sdreport = FALSE,
                              N = 1000,
                              verbose = TRUE,
                              ...) {

  
  # If model is not specified, allow function to choose based on dat_tmb
  # This also abstracts esoteric model specification from the user
  if (is.null(mod)) { # should be it's own function! Can test more easily then
    
    if (!is.null(parameters)) {
      param_names <- names(parameters)
    } else if (!is.null(fit)) {
      param_names <- names(fit$par)
      # add mapped parameters which won't be in fit$par, if appropriate
      if (!is.null(maps)) param_names <- c(param_names, names(maps))
    } else {
      stop("Please provide one of `parameters` or `fit`")
    }
    
    # Start with model with no type information
    mod <- "Surv_SpaceAgeTime"
    
    # determine whether there was information on circumcision type in `out`
    cond <- "type_info" %in% names(dat_tmb) && dat_tmb$type_info == TRUE ||
      # if we don't define dat_tmb, and instead have fit to be re-sampled from
      (!is.null(fit) && any(grepl("mmc", param_names)))
    
    if (cond) {
      mod <- paste0(mod, "_ByType_withUnknownType")
    }
      
    if (any(grepl("paed", param_names))) {
      mod <- paste0(mod, "_Const_Paed_MMC")
    }
    
    # if there are no MMC time autocorr hyperpars, use RW model
    cond <- grepl("logitrho_mmc", param_names) & grepl("time", param_names)
    # if there is a TMC autocorr hyperpar, use RW model for MMC and AR for TMC
    tmc_cond <- "logitrho_tmc_time1" %in% param_names
    if (all(cond == FALSE)) {
      if (tmc_cond) mod <- paste0(mod, "_RW_MMC") else mod <- paste0(mod, "RW")
    }

    # if there is a time term for TMC, use the model with non-constant TMC
    # if (dat_tmb$type_info == TRUE && "u_time_tmc" %in% param_names) {
    cond <- dat_tmb$type_info == TRUE
    if (length(cond) == 0) cond <- "u_time_tmc" %in% param_names
    cond <- cond && "u_time_tmc" %in% param_names
    if (cond) {
      mod <- paste0(mod, "2")
    }
    
    # "RW_MMC" mod only valid where time TMC effect is used
    if (grepl("_RW_MMC", mod) && !grepl("2", mod)) {
      stop(paste(
        "Model with RW for MMC temporal parameters only available when ",
        "including time TMC effect"
      ))
    }
    message("mod not supplied, mod used = ", mod)
  }
  
  if (is.null(mod)) stop("Please provide one of `mod`, `parameters` or `fit`")
  
  # for specified "smaller fit" object (i.e. fit which requires resampling)
  if (!is.null(fit)) {
    if (!is.null(fit$sample)) {
      message("Sample already present in fit object, returning `fit`")
      return(fit)
    }
    if (!is.null(dat_tmb)) {
      message(paste0(
        "No need to specify dat_tmb or parameters for non-null fit, as they",
        " are replaced by those stored in fit"
      ))
    }
    # pull dat_tmb and parameters from small fit
    dat_tmb <- fit$tmb_data
    parameters <- split(fit$par.full, names(fit$par.full))
    init_params <- fit$par_init
    # pull different pars depending on whether the model has mmc/tmc split
    if (!mod %in% c("Surv_SpaceAgeTime", "Surv_SpaceAgeTime_RW")) {
      fit$par_init <- fit$par_init[names(fit$par_init) %in% names(parameters)]
      parameters <- parameters[names(fit$par_init)]
      
      if (any(is.na(names(parameters)))) {
        message("Removing NA parameters, may want to check specifications...")
        parameters <- parameters[!is.na(names(parameters))]
      }
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
        names(init_params) %chin% names(maps)
      ]
      names(parameters)[mapped_pars] <- names(maps)
    }

    is_matrix <- vapply(init_params, is.matrix, logical(1))
    parameters[is_matrix] <- Map(matrix,
      parameters[is_matrix],
      nrow = lapply(init_params[is_matrix], nrow),
      ncol = lapply(init_params[is_matrix], ncol)
    )
    # if no fit == NULL, must have non-null dat_tmb & parameters
  } else {
    if (is.null(dat_tmb) || is.null(parameters)) {
      stop("Please specify non-null dat_tmb and parameters")
    }
  }

  # remove "mmc" from parameter & matrix names if required
  if (mod %in% c("Surv_SpaceAgeTime", "Surv_SpaceAgeTime_RW")) {
    remove_type_distinction <- function(x) {
      names(x) <- stringr::str_remove(names(x), "_mmc")
      x <- x[!grepl("_tmc", names(x))]
    }

    dat_tmb <- remove_type_distinction(
      dat_tmb[!names(dat_tmb) %chin% c("A_mmc", "A_tmc")]
    )
    names(dat_tmb)[names(dat_tmb) == "A_mc"] <- "A"

    parameters <- remove_type_distinction(parameters)
    randoms <- unique(stringr::str_remove(randoms, "_tmc|_mmc"))
  }

  # Only have named random parameters
  randoms <- randoms[randoms %chin% names(parameters)]
  if (length(randoms) == 0) randoms <- NULL
  
  # remove null parameters
  null_pars <- vapply(parameters, is.null, FUN.VALUE = logical(1))
  if (any(null_pars)) {
    message("Removing NULL parameters, check specification...")
    parameters <- parameters[!null_pars]
  }

  if (verbose) message("Creating TMB object with `TMB::MakeADFun`...")
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

    if (verbose) message("Resampling from `fit`...")
    fit$obj <- obj
    fit$obj$fn()
    fit <- circ_sample_tmb(
      fit = fit, obj = obj, nsample = N, sdreport = sdreport
    )
    fit$tmb_data <- fit$par_init <- NULL # make fit object smaller for saving
    return(fit)
  }

  # Run optimiser (use optim if all pars are fixed, nlminb otherwise)
  if (verbose) message("Optimising...")
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
    if (verbose) message("Sampling...")
    fit <- circ_sample_tmb(
      obj = obj, opt = opt, nsample = N, sdreport = sdreport
    )
    # return smaller fit object
    if (smaller_fit_obj == TRUE) {
      if (verbose) message("Minimising fit object...")
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
circ_sample_tmb <- function(fit = NULL,
                            obj = NULL,
                            opt,
                            sdreport = FALSE,
                            nsample = 1000,
                            ...) {

  # Get TMB into "Naomi" format to sample from using naomi
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

  # Generate samples
  fit <- naomi::sample_tmb(fit, nsample = nsample, ...)

  # ensure names for MC columns in fit have the suffix "_mc"
  fit$sample <- append_mc_name(fit$sample)

  return(fit)
}

#### minimise_fit_obj ####

#' @title Minimise Fit Object Size
#' @description Return minimised fit object. Often useful when saving the fit
#' object for later aggregation.
#' @inheritParams threemc_fit_model
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

#### threemc_initial_pars ####

#' @title Initialise `thremec` (hyper)parameters.
#' @description Return minimised fit object. Often useful when saving the fit
#' object for later aggregation.
#' @inheritParams threemc_fit_model
#' @inheritParams threemc_prepare_model_data
#' @param rw_order_tmc_ar Whether to use an AR 1 temporal prior for TMC, 
#' regardless of whether you are using a RW temporal prior for TMC or not, 
#' Default: FALSE
#' @param custom_init named \code{list} of custom fixed and random
#' model parameters you want to supersede "hardcoded" defaults, default = NULL.
#' @return Named \code{list} of intial (hyper)parameters for
#' `threemc_fit_model`
#' @rdname threemc_initial_pars
#' @export
threemc_initial_pars <- function(dat_tmb,
                                 custom_init = NULL,
                                 rw_order = NULL,
                                 rw_order_tmc_ar = FALSE,
                                 paed_age_cutoff = NULL,
                                 inc_time_tmc = FALSE) {


  # Create dummy matrices if not in dat_tmb for particular model specification:

  # dummy paediatric MMC matrices
  if (is.null(paed_age_cutoff)) {
    if ("X_fixed_mmc_paed" %in% names(dat_tmb)) {
      stop(
        "paed_age_cutoff = NULL but dat_tmb$X_fixed_mmc_paed exists, ",
        " please recheck threemc_prepare_model_data arguments"
      )
    }
    X_fixed_mmc_paed <- X_age_mmc_paed <- X_space_mmc_paed <- data.frame(0)
  }

  # dummy time TMC matrices
  if (inc_time_tmc == FALSE) {
    if ("X_time_tmc" %in% names(dat_tmb)) {
      stop(
        "inc_time_tmc = FALSE but dat_tmb$X_time_tmc exists, ",
        " please recheck threemc_prepare_model_data arguments"
      )
    }
    X_time_tmc <- data.frame(0)
  }

  # Initial values
  parameters <- with(
    dat_tmb,
    list(
      # intercept
      "u_fixed_mmc"            = rep(-5, ncol(X_fixed_mmc)),
      "u_fixed_mmc_paed"       = rep(-5, ncol(X_fixed_mmc_paed)),
      "u_fixed_tmc"            = rep(-5, ncol(X_fixed_tmc)),
      # age random effect
      "u_age_mmc"              = rep(0, ncol(X_age_mmc)),
      "u_age_mmc_paed"         = rep(0, ncol(X_age_mmc_paed)),
      "u_age_tmc"              = rep(0, ncol(X_age_tmc)),
      # time random effect for (non-paed) MMC
      "u_time_mmc"             = rep(0, ncol(X_time_mmc)),
      # time random effect for TMC
      "u_time_tmc"             = rep(0, ncol(X_time_tmc)),
      # Space random effect (district)
      "u_space_mmc"            = rep(0, ncol(X_space_mmc)),
      "u_space_mmc_paed"       = rep(0, ncol(X_space_mmc_paed)),
      "u_space_tmc"            = rep(0, ncol(X_space_tmc)),
      # Interactions for MMC
      "u_agetime_mmc"          = matrix(0, ncol(X_age_mmc), ncol(X_time_mmc)),
      "u_agespace_mmc"         = matrix(0, ncol(X_age_mmc), ncol(X_space_mmc)),
      "u_spacetime_mmc"        = matrix(
        0, ncol(X_time_mmc), ncol(X_space_mmc)
      ),
      "u_agespace_mmc_paed"    = matrix(
        0, ncol(X_age_mmc_paed), ncol(X_space_mmc_paed)
      ),
      # Interactions for TMC
      "u_agespace_tmc"         = matrix(0, ncol(X_age_tmc), ncol(X_space_tmc)),
      # Autocorrelation parameters for priors
      # Variance
      "logsigma_age_mmc"            = 0,
      "logsigma_age_mmc_paed"       = 0,
      "logsigma_time_mmc"           = 0,
      "logsigma_space_mmc"          = 0,
      "logsigma_space_mmc_paed"     = 0,
      "logsigma_agetime_mmc"        = 0,
      "logsigma_agespace_mmc"       = 0,
      "logsigma_agespace_mmc_paed"  = 0,
      "logsigma_spacetime_mmc"      = 0,
      "logsigma_age_tmc"            = 0,
      "logsigma_time_tmc"           = 0,
      "logsigma_space_tmc"          = 0,
      "logsigma_agespace_tmc"       = 0,
      # Mean
      "logitrho_mmc_time1"          = 2,
      "logitrho_mmc_time2"          = 2,
      "logitrho_mmc_time3"          = 2,
      "logitrho_mmc_age1"           = 2,
      "logitrho_mmc_paed_age1"      = 2,
      "logitrho_mmc_age2"           = 2,
      "logitrho_mmc_paed_age2"      = 2,
      "logitrho_mmc_age3"           = 2,
      "logitrho_tmc_time1"          = 2,
      "logitrho_tmc_age1"           = 2,
      "logitrho_tmc_age2"           = 2
    )
  )

  # remove paed-related parameters if not desired
  if (is.null(paed_age_cutoff)) {
    parameters <- parameters[!grepl("paed", names(parameters))]
  }

  # remove mmc time correlation parameters, if fitting with RW precision matrix
  # if ("Q_time" %in% names(dat_tmb)) {
  if (!is.null(rw_order)) {
    if (!"Q_time" %in% names(dat_tmb)) {
      stop(
        "rw_order != NULL but dat_tmb$Q_time does not exist, ",
        " please recheck threemc_prepare_model_data arguments"
      )
    }
    parameters <- parameters[!grepl("logitrho_mmc_time", names(parameters))]
    
    # if using an AR 1 temporal prior for TMC, only remove MMC autocorr pars
    remove_pars <- c("logitrho_mmc_time")
    if (rw_order_tmc_ar == FALSE) {
      remove_pars <- c(remove_pars, "logitrho_tmc_time")
    }
    
    parameters <- parameters[!grepl(
      paste(remove_pars, collapse = "|"),
      names(parameters)
    )]
  }

  # remove time tmc terms, if not fitting model with non-constant tmc over time
  if (inc_time_tmc == FALSE) {
    parameters <- parameters[
      !names(parameters) %in% c(
        "u_time_tmc", "logsigma_time_tmc", "logitrho_tmc_time1"
      )
    ]
  }

  # Allow for any custom changes to parameter values
  if (!is.null(custom_init)) {
    parameters[names(parameters) == names(custom_init)] <- custom_init
  }

  return(parameters)
}
