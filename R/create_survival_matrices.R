# function to create survival matrices
create_survival_matrices <- function(
  out, 
  time1 = "time1",
  time2 = "time2",
  age = "age",
  strat = "space") {
  
  out$time1 <- out$time - out$circ_age
  out$time2 <- out$time
  
  #' calculate emppirical agetime hazard matrices for different circ types
  circs <- c(
    "obs_mmc", # medical circumcision rate
    "obs_tmc", # traditional circumcision rate
    "cens",    # censored
    "icens"    # left censored
  )
  hazard_matrices <- lapply(circs, function(x) {
    create_hazard_matrix_agetime(
      dat = out,
      time1 = time1,
      time2 = time2,
      strat = strat,
      age   = age,
      circ  = x,
      Ntime = length(unique(out$time))
    )
  })
  
  # Matrix for selecting instantaneous hazard rate for MMC rate
  A_mmc <- hazard_matrices[[1]]
  # Matrix for selecting instantaneous hazard rate for TMC rate
  A_tmc <- hazard_matrices[[2]]
  # Matrix for selecting instantaneous hazard rate (censored) 
  B <- hazard_matrices[[3]]
  # Matrix for selecting instantaneous hazard rate (left censored)
  C <- hazard_matrices[[4]]
  
  output <- list(
    "A_mmc" = A_mmc,
    "A_tmc" = A_tmc,
    "B"     = B,
    "C"     = C
  )
  return(output)
}