/// @file threemc_no_type.h

/*******************************************************************/
/* Function to calculate nll where we have NO type information     */
/*******************************************************************/
template <class Type>
Type threemc_no_type(
    // sparse model matrices
    const density::SparseMatrix<Type> A,
    const density::SparseMatrix<Type> B, 
    const density::SparseMatrix<Type> C, 
    const density::SparseMatrix<Type> IntMat1, 
    const density::SparseMatrix<Type> IntMat2, 
    
    const density::SparseMatrix<Type> X_fixed, 
    const density::SparseMatrix<Type> X_time,
    const density::SparseMatrix<Type> X_age, 
    const density::SparseMatrix<Type> X_space,
    const density::SparseMatrix<Type> X_agetime, 
    const density::SparseMatrix<Type> X_agespace,
    const density::SparseMatrix<Type> X_spacetime, 
    
    const density::SparseMatrix<Type> Q_space,
    
    vector<Type> u_fixed, vector<Type> u_age,
    vector<Type> u_time, vector<Type> u_space, 
    
    array<Type> u_agetime, array<Type> u_agespace,
    array<Type> u_spacetime, 
    
    Type logsigma_age, Type logsigma_time, Type logsigma_space,
    Type logsigma_agetime, Type logsigma_agespace, Type logsigma_spacetime,
    
    Type logitrho_time1, Type logitrho_time2, Type logitrho_time3,
    Type logitrho_age1, Type logitrho_age2, Type logitrho_age3, 
    
    // reference to struct with report values; function can only return nll
    // report_values<Type>& report_vals
    vector<Type>& haz,    vector<Type>& inc,
    vector<Type>& cum_inc, vector<Type>& surv

){
  
  // Standard deviations
  Type sigma_age       = exp(logsigma_age);
  Type sigma_time      = exp(logsigma_time);
  Type sigma_space     = exp(logsigma_space);
  Type sigma_agetime   = exp(logsigma_agetime);
  Type sigma_agespace  = exp(logsigma_agespace);
  Type sigma_spacetime = exp(logsigma_spacetime);
  
  // Autocorrelation parameters
  Type rho_time1  = geninvlogit(logitrho_time1, Type(-1.0), Type(1.0));
  Type rho_time2  = geninvlogit(logitrho_time2, Type(-1.0), Type(1.0));
  Type rho_time3  = geninvlogit(logitrho_time3, Type(-1.0), Type(1.0));
  Type rho_age1   = geninvlogit(logitrho_age1,  Type(-1.0), Type(1.0));
  Type rho_age2   = geninvlogit(logitrho_age2,  Type(-1.0), Type(1.0));
  Type rho_age3   = geninvlogit(logitrho_age3,  Type(-1.0), Type(1.0));
  
  //////////////////////////////////
  /// Prior on the fixed effects ///
  //////////////////////////////////
  // Negaive log likelihood definition
  Type nll = Type(0);
  
  // Fixed effects for the circumcision rate
  nll -= sum(dnorm(u_fixed,  Type(0), Type(5), TRUE));
  
  ////////////////////////////////////////////
  /// Prior on the temporal random effects ///
  ////////////////////////////////////////////
  // AR1 Process
  nll += density::AR1(rho_time1)(u_time);
  
  // Sum to zero constraint
  nll -= dnorm(u_time.sum(), Type(0), Type(0.001) * u_time.size(), TRUE);
  
  // Prior on the standard deviation for the temporal random effects
  nll -= dexp(sigma_time, Type(1), TRUE) + logsigma_time;
  
  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_time1, Type(3), Type(1), TRUE);
  
  ///////////////////////////////////////
  /// Prior on the age random effects ///
  ///////////////////////////////////////
  // AR1 Process
  nll += density::AR1(rho_age1)(u_age);
  
  // Sum to zero constraint
  nll -= dnorm(u_age.sum(), Type(0), Type(0.001) * u_age.size(), TRUE);
  
  // Prior on the standard deviation for the age random effects
  nll -= dexp(sigma_age, Type(1), TRUE) + logsigma_age;
  
  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_age1, Type(3), Type(1), TRUE);
  
  //////////////////////////////////////////////////
  /// Prior on the stratification random effects ///
  //////////////////////////////////////////////////
  // Gaussian markov random field with prespecified precision matrix
  nll += density::GMRF(Q_space)(u_space);
  
  // Sum to zero constraint
  nll -= dnorm(u_space.sum(), Type(0), Type(0.001) * u_space.size(), TRUE);
  
  // Prior on the standard deviation for the stratification random effects
  nll -= dexp(sigma_space, Type(1), TRUE) + logsigma_space;
  
  ///////////////////////////////////////////////
  /// Prior on the interaction random effects ///
  ///////////////////////////////////////////////
  // Interactions: space-time (GMRF x AR1), age-time (AR1 x AR1) and age-space (AR1 x GMRF)
  nll += SEPARABLE(density::AR1(rho_time2), density::AR1(rho_age2))(u_agetime);
  nll += SEPARABLE(density::GMRF(Q_space), density::AR1(rho_age3))(u_agespace);
  nll += SEPARABLE(density::GMRF(Q_space), density::AR1(rho_time3))(u_spacetime);
  
  // Sum-to-zero constraints
  for (int i = 0; i < u_agetime.cols(); i++) {
    nll -= dnorm(u_agetime.col(i).sum(), Type(0), Type(0.001) * u_agetime.col(i).size(), true);
  }
  for (int i = 0; i < u_agespace.cols(); i++) {
    nll -= dnorm(u_agespace.col(i).sum(), Type(0), Type(0.001) * u_agespace.col(i).size(), true);
  }
  for (int i = 0; i < u_spacetime.cols(); i++) {
    nll -= dnorm(u_spacetime.col(i).sum(), Type(0), Type(0.001) * u_spacetime.col(i).size(), true);
  }
  
  // Vectorising the interaction
  vector<Type> u_agetime_v(u_agetime);
  vector<Type> u_agespace_v(u_agespace);
  vector<Type> u_spacetime_v(u_spacetime);
  
  // Prior on the standard deviations
  nll -= dexp(sigma_agetime, Type(1), TRUE) + logsigma_agetime;
  nll -= dexp(sigma_agespace, Type(1), TRUE) + logsigma_agespace;
  nll -= dexp(sigma_spacetime, Type(1), TRUE) + logsigma_spacetime;
  
  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_time2, Type(3), Type(1), TRUE);
  nll -= dnorm(logitrho_age2, Type(3), Type(1), TRUE);
  nll -= dnorm(logitrho_time3, Type(3), Type(1), TRUE);
  nll -= dnorm(logitrho_age3, Type(3), Type(1), TRUE);
  
  //////////////////////////////
  /// Estimating hazard rate ///
  //////////////////////////////
  // Hazard Rate
 haz = X_fixed * u_fixed +
    X_age * u_age * sigma_age +
    X_space * u_space * sigma_space +
    X_time * u_time * sigma_time +
    X_agetime * u_agetime_v * sigma_agetime +
    X_agespace * u_agespace_v * sigma_agespace +
    X_spacetime * u_spacetime_v * sigma_spacetime;
  
  // Rates on [0,1] scale
  haz = invlogit_vec(haz);
  
  // Survival probabilities
  vector<Type> logprob  = log(Type(1.0) - haz);
               surv     = exp(IntMat1 * logprob);
  vector<Type> surv_lag = exp(IntMat2 * logprob);
  vector<Type> leftcens = Type(1.0) - surv;
  
  // Incidence
  inc = haz * surv_lag;
  
  // Cumulative incidence
  cum_inc = IntMat1 * inc;
  
  //////////////////
  /// Likelihood ///
  //////////////////
  // Getting likelihood for those circumcised
  nll -= (A * log(inc)).sum();
  
  // Getting likelihood for those right censored
  nll -= (B * log(surv)).sum();
  
  // Getting likelihood for those left censored
  nll -= (C * log(leftcens)).sum();
  
  /////////////////////////////////////////
  /// Returning negative log likelihood ///
  /////////////////////////////////////////
  return nll; 
}
