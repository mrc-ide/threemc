
#' @title Produce TMB model fit with sample
#' @description Optimises threemc objective function and produces samples from
#' model fit.
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
#' @param sdreport If set to TRUE, produces the standard deviation report for
#' the model, Default: FALSE
#' @param N Number of samples to be generated, Default: 1000
#' @param ... Further arguments passed to internal functions.
#' @return TMB model fit, including optimised parameters, hessian, samples and
#' standard deviation report (if desired).
#' @rdname threemc_fit_model
#' @export
threemc_fit_model <- function(
  dat_tmb, parameters, maps = NULL,
  randoms = c(
    "u_time_mmc", "u_age_mmc", "u_space_mmc",
    "u_agetime_mmc", "u_agespace_mmc", "u_spacetime_mmc",
    "u_age_tmc", "u_space_tmc", "u_agespace_tmc"
  ),
  mod, sample = TRUE, sdreport = FALSE, N = 1000, ...
) {

  # remove "mmc" from parameter & matrix names if required
  if (mod == "Surv_SpaceAgeTime") {
    remove_type_distinction <- function(x) {
      names(x) <- stringr::str_remove(names(x), "_mmc")
      x <- x[!names(x) %like% "_tmc"]
    }

    dat_tmb <- remove_type_distinction(
      dat_tmb[!names(dat_tmb) %in% c("A_mmc", "A_tmc")]
    )
    names(dat_tmb)[names(dat_tmb) == "A_mc"] <- "A"

    parameters <- remove_type_distinction(parameters)

    randoms <- stringr::str_remove(randoms, "_mmc")
    randoms <- randoms[!randoms %like% "_tmc"]
  }

  # Only have named random parameters
  randoms <- randoms[randoms %in% names(parameters)]
  if (length(randoms) == 0) {
    randoms <- NULL
  }

  # Create TMB object
  obj <- TMB::MakeADFun(dat_tmb,
    parameters,
    random = randoms,
    map = maps,
    method = "BFGS",
    hessian = TRUE,
    DLL = mod,
    ...
  )

  # Run optimiser (use optim if all pars are fixed, nlminb otherwise)
  if (length(maps) == length(obj$par)) {
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
    return(circ_sample_tmb(obj, opt, nsample = N, sdreport = sdreport))
  } else {
    return(opt)
  }
}
