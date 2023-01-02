/// @file threemc_type.h

/************************************************************************/
/* Objective function to specify model and to optimize model parameters */
/* This model splits medical and traditional circumcision types         */
/************************************************************************/
template <class Type>
Type threemc(
    // sparse model matrices
    const density::SparseMatrix<Type> A_mmc,
    const density::SparseMatrix<Type> A_tmc, 
    const density::SparseMatrix<Type> A_mc, 
    const density::SparseMatrix<Type> B, 
    const density::SparseMatrix<Type> C, 
    const density::SparseMatrix<Type> IntMat1, 
    const density::SparseMatrix<Type> IntMat2, 
    
    const density::SparseMatrix<Type> X_fixed_mmc, 
    const density::SparseMatrix<Type> X_time_mmc,
    const density::SparseMatrix<Type> X_age_mmc, 
    const density::SparseMatrix<Type> X_space_mmc,
    const density::SparseMatrix<Type> X_agetime_mmc, 
    const density::SparseMatrix<Type> X_agespace_mmc,
    const density::SparseMatrix<Type> X_spacetime_mmc, 
    // TMC terms
    const density::SparseMatrix<Type> X_fixed_tmc,
    const density::SparseMatrix<Type> X_age_tmc, 
    const density::SparseMatrix<Type> X_space_tmc,
    const density::SparseMatrix<Type> X_agespace_tmc,

    const density::SparseMatrix<Type> Q_space,
  
    vector<Type> u_fixed_mmc, 
    vector<Type> u_fixed_tmc, vector<Type> u_age_mmc,
    vector<Type> u_age_tmc, vector<Type> u_time_mmc, vector<Type> u_space_mmc, 
    vector<Type> u_space_tmc,
    
    array<Type> u_agetime_mmc, array<Type> u_agespace_mmc,
    array<Type> u_spacetime_mmc, array<Type> u_agespace_tmc,
    
    Type logsigma_age_mmc, Type logsigma_time_mmc, Type logsigma_space_mmc,
    Type logsigma_agetime_mmc, Type logsigma_agespace_mmc, Type logsigma_spacetime_mmc,
    Type logsigma_age_tmc, Type logsigma_space_tmc, Type logsigma_agespace_tmc,

    Type logitrho_mmc_time1, Type logitrho_mmc_time2, Type logitrho_mmc_time3,
    Type logitrho_mmc_age1, 
    Type logitrho_mmc_age2, 
    Type logitrho_mmc_age3, Type logitrho_tmc_age1, Type logitrho_tmc_age2,

    // reference to struct with report values; function can only return nll
    // report_values<Type>& report_vals
    vector<Type>& haz_mmc,     vector<Type>& haz_tmc,     vector<Type>& haz, 
    vector<Type>& inc_mmc,     vector<Type>& inc_tmc,     vector<Type>& inc,
    vector<Type>& cum_inc_mmc, vector<Type>& cum_inc_tmc, vector<Type>& cum_inc,
    vector<Type>& surv
  ) {
  
  Type sigma_age_mmc       = exp(logsigma_age_mmc);
  Type sigma_time_mmc      = exp(logsigma_time_mmc);
  Type sigma_space_mmc     = exp(logsigma_space_mmc);
  Type sigma_agetime_mmc   = exp(logsigma_agetime_mmc);
  Type sigma_agespace_mmc  = exp(logsigma_agespace_mmc);
  Type sigma_spacetime_mmc = exp(logsigma_spacetime_mmc);
  Type sigma_age_tmc       = exp(logsigma_age_tmc);
  Type sigma_space_tmc     = exp(logsigma_space_tmc);
  Type sigma_agespace_tmc  = exp(logsigma_agespace_tmc);

  Type rho_mmc_time1  = geninvlogit(logitrho_mmc_time1, Type(-1.0), Type(1.0));
  Type rho_mmc_time2  = geninvlogit(logitrho_mmc_time2, Type(-1.0), Type(1.0));
  Type rho_mmc_time3  = geninvlogit(logitrho_mmc_time3, Type(-1.0), Type(1.0));
  Type rho_mmc_age1   = geninvlogit(logitrho_mmc_age1,  Type(-1.0), Type(1.0));
  Type rho_mmc_age2   = geninvlogit(logitrho_mmc_age2,  Type(-1.0), Type(1.0));
  Type rho_mmc_age3   = geninvlogit(logitrho_mmc_age3,  Type(-1.0), Type(1.0));
  Type rho_tmc_age1   = geninvlogit(logitrho_tmc_age1,  Type(-1.0), Type(1.0));
  Type rho_tmc_age2   = geninvlogit(logitrho_tmc_age2,  Type(-1.0), Type(1.0));

  //////////////////////
  /// Model Variants ///
  /////////////////////
  
  // int rw_order = 0;        // 1 if using RW temporal prior, 0 if using AR 1
  // int paed_age_cutoff = 0; // 1 if including fixed paedatric MMC rate
  // int inc_time_tmc = 0;    // 1 if including time effect for TMC


  //////////////////////////////////
  /// Prior on the fixed effects (make into function!) ///
  //////////////////////////////////  
  
  // Negative log likelihood definition
  Type nll = Type(0);
  
  // Fixed effects for the medical circumcision rate
  nll -= dnorm(u_fixed_mmc,  Type(0), Type(5), TRUE).sum();

  // Fixed effects for the traditional circumcision rate
  nll -= dnorm(u_fixed_tmc, Type(0), Type(5), TRUE).sum();

  // ////////////////////////////////////////////
  // /// Prior on the temporal random effects ///
  // ////////////////////////////////////////////
  // AR1 Process
  nll += density::AR1(rho_mmc_time1)(u_time_mmc);

  // Sum to zero constraint
  nll -= dnorm(u_time_mmc.sum(), Type(0), Type(0.001) * u_time_mmc.size(), TRUE);

  // Prior on the standard deviation for the temporal random effects
  nll -= dexp(sigma_time_mmc, Type(1), TRUE) + logsigma_time_mmc;

  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_mmc_time1, Type(3), Type(3), TRUE);

  ///////////////////////////////////////
  /// Prior on the age random effects ///
  ///////////////////////////////////////
  // density::AR1 processes
  nll += density::AR1(rho_mmc_age1)(u_age_mmc);
  nll += density::AR1(rho_tmc_age1)(u_age_tmc);

  // Sum to zero constraint
  nll -= dnorm(u_age_mmc.sum(), Type(0), Type(0.001) * u_age_mmc.size(), TRUE);
  nll -= dnorm(u_age_tmc.sum(), Type(0), Type(0.001) * u_age_tmc.size(), TRUE);

  // Prior on the standard deviation for the age random effects
  nll -= dexp(sigma_age_mmc, Type(1), TRUE) + logsigma_age_mmc;
  nll -= dexp(sigma_age_tmc, Type(1), TRUE) + logsigma_age_tmc;

  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_mmc_age1, Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_tmc_age1, Type(3), Type(2), TRUE);

  ///////////////////////////////////////////
  /// Prior on the spatial random effects ///
  ///////////////////////////////////////////
  // Gaussian markov random field with prespecified precision matrix
  nll += density::GMRF(Q_space)(u_space_mmc);
  nll += density::GMRF(Q_space)(u_space_tmc);

  // Sum to zero constraints
  nll -= dnorm(u_space_mmc.sum(), Type(0), Type(0.001) * u_space_mmc.size(), TRUE);
  nll -= dnorm(u_space_tmc.sum(), Type(0), Type(0.001) * u_space_tmc.size(), TRUE);

  // Prior on the standard deviation for the spatial random effects
  nll -= dexp(sigma_space_mmc, Type(1), TRUE) + logsigma_space_mmc;
  nll -= dexp(sigma_space_tmc, Type(1), TRUE) + logsigma_space_tmc;

  ///////////////////////////////////////////////
  /// Prior on the interaction random effects ///
  ///////////////////////////////////////////////
  // Interactions: space-time (GMRF x AR1), age-time (AR1 x AR1) and age-space (AR1 x GMRF)
  nll += SEPARABLE(density::AR1(rho_mmc_time2), density::AR1(rho_mmc_age2))(u_agetime_mmc);
  nll += SEPARABLE(density::GMRF(Q_space), density::AR1(rho_mmc_age3))(u_agespace_mmc);
  nll += SEPARABLE(density::GMRF(Q_space), density::AR1(rho_mmc_time3))(u_spacetime_mmc);
  nll += SEPARABLE(density::GMRF(Q_space), density::AR1(rho_tmc_age2))(u_agespace_tmc);
  
  // Sum-to-zero constraints
  for (int i = 0; i < u_agespace_mmc.cols(); i++) {
    nll -= dnorm(u_agespace_mmc.col(i).sum(), Type(0), Type(0.001) * u_agespace_mmc.col(i).size(), true);
  } 
  for (int i = 0; i < u_agetime_mmc.cols(); i++) {
    nll -= dnorm(u_agetime_mmc.col(i).sum(), Type(0), Type(0.001) * u_agetime_mmc.col(i).size(), true);
  }  
  for (int i = 0; i < u_spacetime_mmc.cols(); i++) {
    nll -= dnorm(u_spacetime_mmc.col(i).sum(), Type(0), Type(0.001) * u_spacetime_mmc.col(i).size(), true);
  }  
  for (int i = 0; i < u_agespace_tmc.cols(); i++) {
    nll -= dnorm(u_agespace_tmc.col(i).sum(), Type(0), Type(0.001) * u_agespace_tmc.col(i).size(), true);
  }  
  
  // Vectorising the interaction
  vector<Type> u_agespace_mmc_v(u_agespace_mmc);
  vector<Type> u_agetime_mmc_v(u_agetime_mmc);
  vector<Type> u_spacetime_mmc_v(u_spacetime_mmc);
  vector<Type> u_agespace_tmc_v(u_agespace_tmc);

  // Prior on the standard deviation for the interaction random effects
  nll -= dexp(sigma_agespace_mmc,  Type(1), TRUE) + logsigma_agespace_mmc;
  nll -= dexp(sigma_agetime_mmc,   Type(1), TRUE) + logsigma_agetime_mmc;
  nll -= dexp(sigma_spacetime_mmc, Type(1), TRUE) + logsigma_spacetime_mmc;
  nll -= dexp(sigma_agespace_tmc,  Type(1), TRUE) + logsigma_agespace_tmc;

  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_mmc_time2, Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_mmc_age2,  Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_mmc_time3, Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_mmc_age3,  Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_tmc_age2,  Type(3), Type(2), TRUE);

  //////////////////////////////
  /// Estimating hazard rate ///
  //////////////////////////////
  // Medical hazard rate
  // vector<Type> haz_mmc = X_fixed_mmc * u_fixed_mmc +
  haz_mmc = X_fixed_mmc * u_fixed_mmc +
    X_time_mmc * u_time_mmc * sigma_time_mmc +
    X_space_mmc * u_space_mmc * sigma_space_mmc +
    X_age_mmc * u_age_mmc * sigma_age_mmc +
    X_agetime_mmc * u_agetime_mmc_v * sigma_agetime_mmc +
    X_agespace_mmc * u_agespace_mmc_v * sigma_agespace_mmc +
    X_spacetime_mmc * u_spacetime_mmc_v * sigma_spacetime_mmc;

  // Traditional hazard rate
  haz_tmc = X_fixed_tmc * u_fixed_tmc +
    X_space_tmc * u_space_tmc * sigma_space_tmc +
    X_age_tmc * u_age_tmc * sigma_age_tmc +
    X_agespace_tmc * u_agespace_tmc_v * sigma_agespace_tmc;

  // Rates on [0,1] scale
  haz_tmc = invlogit_vec(haz_tmc);
  haz_mmc = invlogit_vec(haz_mmc);

  // Adjustment such that \lambda_mmc + \lambda_tmc \in [0,1]
  // Medical rate to only take from the remaining proportion
  // not taken through traditional circumcision (1 - \lambda_tmc)
  haz_mmc = haz_mmc * (1 - haz_tmc);

  // Total hazard rate
  haz = haz_mmc + haz_tmc;

  // Survival probabilities
  vector<Type> logprob  = log(Type(1.0) - haz);
  surv     = exp(IntMat1 * logprob);
  vector<Type> surv_lag = exp(IntMat2 * logprob);
  vector<Type> leftcens = Type(1.0) - surv;

  // Incidence
  inc_tmc = haz_tmc * surv_lag;
  inc_mmc = haz_mmc * surv_lag;
  inc = haz * surv_lag;

  // Cumulative incidence
  cum_inc_tmc = IntMat1 * inc_tmc;
  cum_inc_mmc = IntMat1 * inc_mmc;
  cum_inc = cum_inc_tmc + cum_inc_mmc;
  
  //////////////////
  /// Likelihood ///
  //////////////////
  // Getting likelihood for those medically circumcised
  nll -= (A_mmc * log(inc_mmc)).sum();

  // Getting likelihood for those traditionally circumcised
  nll -= (A_tmc * log(inc_tmc)).sum();

  // Getting likelihood for those circumcised of unknown type
  nll -= (A_mc * log(inc)).sum();

  // Getting likelihood for those right censored
  nll -= (B * log(surv)).sum();

  // Getting likelihood for those left censored
  nll -= (C * log(leftcens)).sum();

  /////////////////////////////////////////
  /// Returning negative log likelihood ///
  /////////////////////////////////////////
  return nll;
}

/************************************************************************/
/* Objective function to specify model and to optimize model parameters */
/* This model splits medical and traditional circumcision types, and    */
/* incorporates a paediatric medical circumcision type.                 */
/************************************************************************/
template <class Type>
Type threemc(
    // sparse model matrices
    const density::SparseMatrix<Type> A_mmc,
    const density::SparseMatrix<Type> A_tmc, 
    const density::SparseMatrix<Type> A_mc, 
    const density::SparseMatrix<Type> B, 
    const density::SparseMatrix<Type> C, 
    const density::SparseMatrix<Type> IntMat1, 
    const density::SparseMatrix<Type> IntMat2, 
    
    // MMC terms
    const density::SparseMatrix<Type> X_fixed_mmc, 
    const density::SparseMatrix<Type> X_time_mmc,
    const density::SparseMatrix<Type> X_age_mmc, 
    const density::SparseMatrix<Type> X_space_mmc,
    const density::SparseMatrix<Type> X_agetime_mmc, 
    const density::SparseMatrix<Type> X_agespace_mmc,
    const density::SparseMatrix<Type> X_spacetime_mmc, 
    // paed MMC terms
    const density::SparseMatrix<Type> X_fixed_mmc_paed, 
    const density::SparseMatrix<Type> X_age_mmc_paed,
    const density::SparseMatrix<Type> X_space_mmc_paed,
    const density::SparseMatrix<Type> X_agespace_mmc_paed,
    // TMC terms
    const density::SparseMatrix<Type> X_fixed_tmc,
    const density::SparseMatrix<Type> X_age_tmc, 
    const density::SparseMatrix<Type> X_space_tmc,
    const density::SparseMatrix<Type> X_agespace_tmc,

    const density::SparseMatrix<Type> Q_space,
  
    vector<Type> u_fixed_mmc, 
    vector<Type> u_fixed_mmc_paed, 
    vector<Type> u_fixed_tmc, vector<Type> u_age_mmc,
    vector<Type> u_age_mmc_paed,
    vector<Type> u_age_tmc, vector<Type> u_time_mmc, vector<Type> u_space_mmc, 
    vector<Type> u_space_mmc_paed,
    vector<Type> u_space_tmc,
    
    array<Type> u_agetime_mmc, array<Type> u_agespace_mmc,
    array<Type> u_agespace_mmc_paed,
    array<Type> u_spacetime_mmc, array<Type> u_agespace_tmc,
    
    Type logsigma_age_mmc, Type logsigma_time_mmc, Type logsigma_space_mmc,
    Type logsigma_age_mmc_paed, Type logsigma_space_mmc_paed,
    Type logsigma_agetime_mmc, Type logsigma_agespace_mmc, Type logsigma_spacetime_mmc,
    Type logsigma_agespace_mmc_paed,
    Type logsigma_age_tmc, Type logsigma_space_tmc, Type logsigma_agespace_tmc,

    Type logitrho_mmc_time1, Type logitrho_mmc_time2, Type logitrho_mmc_time3,
    Type logitrho_mmc_age1, Type logitrho_mmc_paed_age1,  
    Type logitrho_mmc_age2, Type logitrho_mmc_paed_age2,    
    Type logitrho_mmc_age3, Type logitrho_tmc_age1, Type logitrho_tmc_age2,

    // reference to struct with report values; function can only return nll
    // report_values<Type>& report_vals
    vector<Type>& haz_mmc,     vector<Type>& haz_tmc,     vector<Type>& haz, 
    vector<Type>& inc_mmc,     vector<Type>& inc_tmc,     vector<Type>& inc,
    vector<Type>& cum_inc_mmc, vector<Type>& cum_inc_tmc, vector<Type>& cum_inc,
    vector<Type>& surv
  ) {
  
  Type sigma_age_mmc       = exp(logsigma_age_mmc);
  Type sigma_age_mmc_paed  = exp(logsigma_age_mmc_paed);
  Type sigma_time_mmc      = exp(logsigma_time_mmc);
  Type sigma_space_mmc     = exp(logsigma_space_mmc);
  Type sigma_space_mmc_paed  = exp(logsigma_space_mmc_paed);
  Type sigma_agetime_mmc   = exp(logsigma_agetime_mmc);
  Type sigma_agespace_mmc  = exp(logsigma_agespace_mmc);
  Type sigma_agespace_mmc_paed = exp(logsigma_agespace_mmc_paed);
  Type sigma_spacetime_mmc = exp(logsigma_spacetime_mmc);
  Type sigma_age_tmc       = exp(logsigma_age_tmc);
  Type sigma_space_tmc     = exp(logsigma_space_tmc);
  Type sigma_agespace_tmc  = exp(logsigma_agespace_tmc);

  Type rho_mmc_time1  = geninvlogit(logitrho_mmc_time1, Type(-1.0), Type(1.0));
  Type rho_mmc_time2  = geninvlogit(logitrho_mmc_time2, Type(-1.0), Type(1.0));
  Type rho_mmc_time3  = geninvlogit(logitrho_mmc_time3, Type(-1.0), Type(1.0));
  Type rho_mmc_age1   = geninvlogit(logitrho_mmc_age1,  Type(-1.0), Type(1.0));
  Type rho_mmc_paed_age1 = geninvlogit(logitrho_mmc_paed_age1,  Type(-1.0), Type(1.0));
  Type rho_mmc_age2   = geninvlogit(logitrho_mmc_age2,  Type(-1.0), Type(1.0));
  Type rho_mmc_paed_age2 = geninvlogit(logitrho_mmc_paed_age2,  Type(-1.0), Type(1.0));
  Type rho_mmc_age3   = geninvlogit(logitrho_mmc_age3,  Type(-1.0), Type(1.0));
  Type rho_tmc_age1   = geninvlogit(logitrho_tmc_age1,  Type(-1.0), Type(1.0));
  Type rho_tmc_age2   = geninvlogit(logitrho_tmc_age2,  Type(-1.0), Type(1.0));

  //////////////////////
  /// Model Variants ///
  /////////////////////
  
  // int rw_order = 0;        // 1 if using RW temporal prior, 0 if using AR 1
  // int paed_age_cutoff = 0; // 1 if including fixed paedatric MMC rate
  // int inc_time_tmc = 0;    // 1 if including time effect for TMC


  //////////////////////////////////
  /// Prior on the fixed effects (make into function!) ///
  //////////////////////////////////  
  
  // Negative log likelihood definition
  Type nll = Type(0);
  
  // Fixed effects for the medical circumcision rate
  nll -= dnorm(u_fixed_mmc,  Type(0), Type(5), TRUE).sum();

  // Fixed effects for the paed medical circumcision rate
  nll -= dnorm(u_fixed_mmc_paed, Type(0), Type(5), TRUE).sum();

  // Fixed effects for the traditional circumcision rate
  nll -= dnorm(u_fixed_tmc, Type(0), Type(5), TRUE).sum();

  // ////////////////////////////////////////////
  // /// Prior on the temporal random effects ///
  // ////////////////////////////////////////////
  // AR1 Process
  nll += density::AR1(rho_mmc_time1)(u_time_mmc);

  // Sum to zero constraint
  nll -= dnorm(u_time_mmc.sum(), Type(0), Type(0.001) * u_time_mmc.size(), TRUE);

  // Prior on the standard deviation for the temporal random effects
  nll -= dexp(sigma_time_mmc, Type(1), TRUE) + logsigma_time_mmc;

  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_mmc_time1, Type(3), Type(3), TRUE);

  ///////////////////////////////////////
  /// Prior on the age random effects ///
  ///////////////////////////////////////
  // density::AR1 processes
  nll += density::AR1(rho_mmc_age1)(u_age_mmc);
  nll += density::AR1(rho_mmc_paed_age1)(u_age_mmc_paed);
  nll += density::AR1(rho_tmc_age1)(u_age_tmc);

  // Sum to zero constraint
  nll -= dnorm(u_age_mmc.sum(), Type(0), Type(0.001) * u_age_mmc.size(), TRUE);
  nll -= dnorm(u_age_mmc_paed.sum(), Type(0), Type(0.001) * u_age_mmc_paed.size(), TRUE);
  nll -= dnorm(u_age_tmc.sum(), Type(0), Type(0.001) * u_age_tmc.size(), TRUE);

  // Prior on the standard deviation for the age random effects
  nll -= dexp(sigma_age_mmc, Type(1), TRUE) + logsigma_age_mmc;
  nll -= dexp(sigma_age_mmc_paed, Type(1), TRUE) + logsigma_age_mmc_paed;
  nll -= dexp(sigma_age_tmc, Type(1), TRUE) + logsigma_age_tmc;

  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_mmc_age1, Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_mmc_paed_age1, Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_tmc_age1, Type(3), Type(2), TRUE);

  ///////////////////////////////////////////
  /// Prior on the spatial random effects ///
  ///////////////////////////////////////////
  // Gaussian markov random field with prespecified precision matrix
  nll += density::GMRF(Q_space)(u_space_mmc);
  nll += density::GMRF(Q_space)(u_space_mmc_paed);
  nll += density::GMRF(Q_space)(u_space_tmc);

  // Sum to zero constraints
  nll -= dnorm(u_space_mmc.sum(), Type(0), Type(0.001) * u_space_mmc.size(), TRUE);
  nll -= dnorm(u_space_mmc_paed.sum(), Type(0), Type(0.001) * u_space_mmc_paed.size(), TRUE);
  nll -= dnorm(u_space_tmc.sum(), Type(0), Type(0.001) * u_space_tmc.size(), TRUE);

  // Prior on the standard deviation for the spatial random effects
  nll -= dexp(sigma_space_mmc, Type(1), TRUE) + logsigma_space_mmc;
  nll -= dexp(sigma_space_mmc_paed, Type(1), TRUE) + logsigma_space_mmc_paed;
  nll -= dexp(sigma_space_tmc, Type(1), TRUE) + logsigma_space_tmc;

  ///////////////////////////////////////////////
  /// Prior on the interaction random effects ///
  ///////////////////////////////////////////////
  // Interactions: space-time (GMRF x AR1), age-time (AR1 x AR1) and age-space (AR1 x GMRF)
  nll += SEPARABLE(density::AR1(rho_mmc_time2), density::AR1(rho_mmc_age2))(u_agetime_mmc);
  nll += SEPARABLE(density::GMRF(Q_space), density::AR1(rho_mmc_age3))(u_agespace_mmc);
  nll += SEPARABLE(density::GMRF(Q_space), density::AR1(rho_mmc_time3))(u_spacetime_mmc);
  nll += SEPARABLE(density::GMRF(Q_space), density::AR1(rho_mmc_paed_age2))(u_agespace_mmc_paed);
  nll += SEPARABLE(density::GMRF(Q_space), density::AR1(rho_tmc_age2))(u_agespace_tmc);
  
  // Sum-to-zero constraints
  for (int i = 0; i < u_agespace_mmc.cols(); i++) {
    nll -= dnorm(u_agespace_mmc.col(i).sum(), Type(0), Type(0.001) * u_agespace_mmc.col(i).size(), true);
  } 
  for (int i = 0; i < u_agetime_mmc.cols(); i++) {
    nll -= dnorm(u_agetime_mmc.col(i).sum(), Type(0), Type(0.001) * u_agetime_mmc.col(i).size(), true);
  }  
  for (int i = 0; i < u_spacetime_mmc.cols(); i++) {
    nll -= dnorm(u_spacetime_mmc.col(i).sum(), Type(0), Type(0.001) * u_spacetime_mmc.col(i).size(), true);
  }  
  for (int i = 0; i < u_agespace_mmc_paed.cols(); i++) {
    nll -= dnorm(u_agespace_mmc_paed.col(i).sum(), Type(0), Type(0.001) * u_agespace_mmc_paed.col(i).size(), true);
  }
  for (int i = 0; i < u_agespace_tmc.cols(); i++) {
    nll -= dnorm(u_agespace_tmc.col(i).sum(), Type(0), Type(0.001) * u_agespace_tmc.col(i).size(), true);
  }  
  
  // Vectorising the interaction
  vector<Type> u_agespace_mmc_v(u_agespace_mmc);
  vector<Type> u_agetime_mmc_v(u_agetime_mmc);
  vector<Type> u_spacetime_mmc_v(u_spacetime_mmc);
  vector<Type> u_agespace_mmc_paed_v(u_agespace_mmc_paed);
  vector<Type> u_agespace_tmc_v(u_agespace_tmc);

  // Prior on the standard deviation for the interaction random effects
  nll -= dexp(sigma_agespace_mmc,  Type(1), TRUE) + logsigma_agespace_mmc;
  nll -= dexp(sigma_agetime_mmc,   Type(1), TRUE) + logsigma_agetime_mmc;
  nll -= dexp(sigma_spacetime_mmc, Type(1), TRUE) + logsigma_spacetime_mmc;
  nll -= dexp(sigma_agespace_mmc_paed, Type(1), TRUE) + logsigma_agespace_mmc_paed;
  nll -= dexp(sigma_agespace_tmc,  Type(1), TRUE) + logsigma_agespace_tmc;

  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_mmc_time2, Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_mmc_age2,  Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_mmc_time3, Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_mmc_age3,  Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_mmc_paed_age2, Type(3), Type(2), TRUE);
  nll -= dnorm(logitrho_tmc_age2,  Type(3), Type(2), TRUE);

  //////////////////////////////
  /// Estimating hazard rate ///
  //////////////////////////////
  // Medical hazard rate
  // vector<Type> haz_mmc = X_fixed_mmc * u_fixed_mmc +
  haz_mmc = X_fixed_mmc * u_fixed_mmc +
    X_time_mmc * u_time_mmc * sigma_time_mmc +
    X_space_mmc * u_space_mmc * sigma_space_mmc +
    X_age_mmc * u_age_mmc * sigma_age_mmc +
    X_agetime_mmc * u_agetime_mmc_v * sigma_agetime_mmc +
    X_agespace_mmc * u_agespace_mmc_v * sigma_agespace_mmc +
    X_spacetime_mmc * u_spacetime_mmc_v * sigma_spacetime_mmc +
    // Peadiatric part of the MMC process 
    X_fixed_mmc_paed * u_fixed_mmc_paed +
    X_space_mmc_paed * u_space_mmc_paed * sigma_space_mmc_paed + 
		X_age_mmc_paed * u_age_mmc_paed * sigma_age_mmc_paed + 
		X_agespace_mmc_paed * u_agespace_mmc_paed_v * sigma_agespace_mmc_paed;

  // Traditional hazard rate
  haz_tmc = X_fixed_tmc * u_fixed_tmc +
    X_space_tmc * u_space_tmc * sigma_space_tmc +
    X_age_tmc * u_age_tmc * sigma_age_tmc +
    X_agespace_tmc * u_agespace_tmc_v * sigma_agespace_tmc;

  // Rates on [0,1] scale
  haz_tmc = invlogit_vec(haz_tmc);
  haz_mmc = invlogit_vec(haz_mmc);

  // Adjustment such that \lambda_mmc + \lambda_tmc \in [0,1]
  // Medical rate to only take from the remaining proportion
  // not taken through traditional circumcision (1 - \lambda_tmc)
  haz_mmc = haz_mmc * (1 - haz_tmc);

  // Total hazard rate
  haz = haz_mmc + haz_tmc;

  // Survival probabilities
  vector<Type> logprob  = log(Type(1.0) - haz);
  surv     = exp(IntMat1 * logprob);
  vector<Type> surv_lag = exp(IntMat2 * logprob);
  vector<Type> leftcens = Type(1.0) - surv;

  // Incidence
  inc_tmc = haz_tmc * surv_lag;
  inc_mmc = haz_mmc * surv_lag;
  inc = haz * surv_lag;

  // Cumulative incidence
  cum_inc_tmc = IntMat1 * inc_tmc;
  cum_inc_mmc = IntMat1 * inc_mmc;
  cum_inc = cum_inc_tmc + cum_inc_mmc;
  
  //////////////////
  /// Likelihood ///
  //////////////////
  // Getting likelihood for those medically circumcised
  nll -= (A_mmc * log(inc_mmc)).sum();

  // Getting likelihood for those traditionally circumcised
  nll -= (A_tmc * log(inc_tmc)).sum();

  // Getting likelihood for those circumcised of unknown type
  nll -= (A_mc * log(inc)).sum();

  // Getting likelihood for those right censored
  nll -= (B * log(surv)).sum();

  // Getting likelihood for those left censored
  nll -= (C * log(leftcens)).sum();

  /////////////////////////////////////////
  /// Returning negative log likelihood ///
  /////////////////////////////////////////
  return nll;
}
