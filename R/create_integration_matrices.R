
# function to create integration matrices for selecting instantaneous hazard
# rate
create_integration_matrices <- function(out,
                                        time1 = "time1",
                                        time2 = "time2",
                                        age = "age",
                                        strat = "space"
                                        ) {
  
  out <- out
  out$time1 <- out$time - out$circ_age
  out$time2 <- out$time
  out$age <- out$circ_age + 1
  
  # Matrix for selecting instantaneous hazard rate
  IntMat1 <- create_integration_matrix_agetime(
    dat = out,
    time1 = time1,
    time2 = time2,
    strat = strat,
    age = age,
    Ntime = length(unique(out$time))
  )
  
  # Matrix for selecting instantaneous hazard rate
  IntMat2 <- create_integration_matrix_agetime_lag(
    dat = out,
    time1 = "time1",
    time2 = "time2",
    strat = "space",
    age = "age",
    Ntime = length(unique(out$time))
  )
  
  output <- list(
    "IntMat1" = IntMat1,
    "IntMat2" = IntMat2
  )
  return(output)
}
