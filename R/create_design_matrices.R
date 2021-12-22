#' @title Create Design Matrices
#' @description Create design matrices for fixed effects and temporal, age, 
#' spatial and random effects, for both medical and traditional circumcision.
#' 
#' @param out Shell dataset (outputted by \link[threemc]{create_shell_dataset}
#' with a row for every unique record in circumcision survey data for a given 
#' area. Also includes empirical estimates for circumcision estimates for each 
#' unique record.
#' @param k_dt Age knot spacing in spline definitions, Default: 5 
#' @return List of design matrices for fixed and random effects for medical 
#' and traditional circumcision.
#
#' @seealso 
#'  \code{\link}[threemc]{create_shell_dataset}}
#'  \code{\link[splines]{splineDesign}}
#'  \code{\link[mgcv]{tensor.prod.model.matrix}}
#'  \code{\link[Matrix]{sparse.model.matrix}}
#' @rdname create_design_matrices
#' @export 
#' 
#' @importFrom splines splineDesign
#' @importFrom mgcv tensor.prod.model.matrix
#' @importFrom Matrix sparse.model.matrix

# function to create design matrices
create_design_matrices <- function(out, k_dt = 5) {
  
  # Spline definitions
  k_age <- 
    k_dt * (floor(min(out$age) / k_dt) - 3):(ceiling(max(out$age) / k_dt) + 3)
  
  # Design matrix for the fixed effects
  X_fixed <- sparse.model.matrix(N ~ 1, data = out)
  
  # Design matrix for the temporal random effects
  X_time <- sparse.model.matrix(N ~ -1 + as.factor(time), data = out)
  
  # Design matrix for the age random effects
  X_age <- splineDesign(k_age, out$age, outer.ok = TRUE)
  X_age <- as(X_age, "sparseMatrix")
  
  # Design matrix for the spatial random effects
  X_space <- sparse.model.matrix(N ~ -1 + as.factor(space), data = out)
  
  # Design matrix for the interaction random effects
  X_agetime <- tensor.prod.model.matrix(list(X_time, X_age))
  X_agespace <- tensor.prod.model.matrix(list(X_space, X_age))
  X_spacetime <- sparse.model.matrix(
    N ~ -1 + factor((out %>% 
                       group_by(space, time) %>%
                       group_indices())), 
    data = out
  )
  # return design matrices as list
  output <- list(
    "X_fixed_mmc"     = X_fixed,
    "X_fixed_tmc"     = X_fixed,
    "X_time_mmc"      = X_time,
    "X_age_mmc"       = X_age,
    "X_age_tmc"       = X_age,
    "X_space_mmc"     = X_space,
    "X_space_tmc"     = X_space,
    "X_agetime_mmc"   = X_agetime,
    "X_agespace_mmc"  = X_agespace,
    "X_agespace_tmc"  = X_agespace,
    "X_spacetime_mmc" = X_spacetime
  )
  return(output)
}
