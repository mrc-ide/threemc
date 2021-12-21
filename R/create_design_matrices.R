# function to create design matrices
create_design_matrices <- function(out, k_dt) {
  
  # Spline definitions
  k_age <- 
    k_dt * (floor(min(out$age) / k_dt) - 3):(ceiling(max(out$age) / k_dt) + 3)
  
  # Design matrix for the fixed effects
  X_fixed <- sparse.model.matrix(N ~ 1, data = out)
  
  # Design matrix for the temporal random effects
  X_time <- sparse.model.matrix(N ~ -1 + as.factor(time), data = out)
  
  # Design matrix for the age random effects
  X_age <- splines::splineDesign(k_age, out$age, outer.ok = TRUE)
  X_age <- as(X_age, "sparseMatrix")
  
  # Design matrix for the spatial random effects
  X_space <- sparse.model.matrix(N ~ -1 + as.factor(space), data = out)
  
  # Design matrix for the interaction random effects
  X_agetime <- mgcv::tensor.prod.model.matrix(list(X_time, X_age))
  X_agespace <- mgcv::tensor.prod.model.matrix(list(X_space, X_age))
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
