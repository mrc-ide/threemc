/// @file threemc.cpp

#define TMB_LIB_INIT R_init_threemc
#include <TMB.hpp>
#include "utils.h"
#include "threemc_type.h"

/************************************************************************/
/* Objective function to specify model and to optimize model parameters */
/************************************************************************/
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  
  using namespace density;
  
  ////////////////////////
  /// Data definitions ///
  ////////////////////////
  // Survival analysis matrices
  DATA_SPARSE_MATRIX(A_mmc); // Matrix selecting instantaneous hazard for medically circumcised pop
  DATA_SPARSE_MATRIX(A_tmc); // Matrix selecting instantaneous hazard for traditionally circumcised pop
  DATA_SPARSE_MATRIX(A_mc); // Matrix selecting instantaneous hazard for unknown circumcised pop
  DATA_SPARSE_MATRIX(B); // Matrix selecting relevant cumulative hazard entry for observed and right censored pop
  DATA_SPARSE_MATRIX(C); // Matrix selecting relevant cumulative hazard entry for interval censored pop
  DATA_SPARSE_MATRIX(IntMat1); // Integration matrix for cumulative hazard 
  DATA_SPARSE_MATRIX(IntMat2); // Integration matrix for lagged cumulative hazard 
  
  // Design matrices 
  DATA_SPARSE_MATRIX(X_fixed_mmc); // Design matrix for the fixed effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_time_mmc); // Design matrix for the temporal random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_age_mmc); // Design matrix for the stratification random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_space_mmc); // Design matrix for the stratification random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_agetime_mmc); // Design matrix for the interaction random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_agespace_mmc); // Design matrix for the interaction random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_spacetime_mmc); // Design matrix for the interaction random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_fixed_tmc); // Design matrix for the fixed effects in the traditional circumcision hazard rate
  DATA_SPARSE_MATRIX(X_age_tmc); // Design matrix for the stratification random effects in the traditional circumcision hazard rate
  DATA_SPARSE_MATRIX(X_space_tmc); // Design matrix for the stratification random effects in the medical circumcision hazard rate
  DATA_SPARSE_MATRIX(X_agespace_tmc); // Design matrix for the interaction random effects in the medical circumcision hazard rate
  
  // Precision matrices 
  DATA_SPARSE_MATRIX(Q_space); // Aggregation matrix for number of circumcisions performed
  
  //////////////////
  /// Parameters ///
  //////////////////
  // Fixed Effects
  PARAMETER_VECTOR(u_fixed_mmc);
  PARAMETER_VECTOR(u_fixed_tmc);
  
  // Age random effect
  PARAMETER_VECTOR(u_age_mmc); 
  PARAMETER_VECTOR(u_age_tmc); 
  
  // Temporal random effects 
  PARAMETER_VECTOR(u_time_mmc);
  
  // Spatial random effects
  PARAMETER_VECTOR(u_space_mmc);
  PARAMETER_VECTOR(u_space_tmc);
  
  // Interactions
  PARAMETER_ARRAY(u_agetime_mmc);
  PARAMETER_ARRAY(u_agespace_mmc);
  PARAMETER_ARRAY(u_spacetime_mmc);
  PARAMETER_ARRAY(u_agespace_tmc);
  
  // Standard deviations 
  PARAMETER(logsigma_age_mmc);       // Type sigma_age_mmc       = exp(logsigma_age_mmc);
  PARAMETER(logsigma_time_mmc);      // Type sigma_time_mmc      = exp(logsigma_time_mmc);
  PARAMETER(logsigma_space_mmc);     // Type sigma_space_mmc     = exp(logsigma_space_mmc);
  PARAMETER(logsigma_agetime_mmc);   // Type sigma_agetime_mmc   = exp(logsigma_agetime_mmc);
  PARAMETER(logsigma_agespace_mmc);  // Type sigma_agespace_mmc  = exp(logsigma_agespace_mmc);
  PARAMETER(logsigma_spacetime_mmc); // Type sigma_spacetime_mmc = exp(logsigma_spacetime_mmc);
  PARAMETER(logsigma_age_tmc);       // Type sigma_age_tmc       = exp(logsigma_age_tmc);
  PARAMETER(logsigma_space_tmc);     // Type sigma_space_tmc     = exp(logsigma_space_tmc);
  PARAMETER(logsigma_agespace_tmc);  // Type sigma_agespace_tmc  = exp(logsigma_agespace_tmc);
  
  // Autocorrelation parameters 
  PARAMETER(logitrho_mmc_time1);  // Type rho_mmc_time1  = geninvlogit(logitrho_mmc_time1, Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_time2);  // Type rho_mmc_time2  = geninvlogit(logitrho_mmc_time2, Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_time3);  // Type rho_mmc_time3  = geninvlogit(logitrho_mmc_time3, Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_age1);   // Type rho_mmc_age1   = geninvlogit(logitrho_mmc_age1,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_age2);   // Type rho_mmc_age2   = geninvlogit(logitrho_mmc_age2,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_age3);   // Type rho_mmc_age3   = geninvlogit(logitrho_mmc_age3,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_tmc_age1);   // Type rho_tmc_age1   = geninvlogit(logitrho_tmc_age1,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_tmc_age2);   // Type rho_tmc_age2   = geninvlogit(logitrho_tmc_age2,  Type(-1.0), Type(1.0));
 
  // Negative log likelihood definition
  Type nll = Type(0);
  
  // Define covariance structure for the conditional model
  DATA_STRUCT(report_vals, report_values);
  
  nll = threemc_type(
    
    A_mmc, A_tmc, A_mc, B, C, IntMat1, IntMat2, 
    
    X_fixed_mmc, X_time_mmc, X_age_mmc, X_space_mmc, X_agetime_mmc,
    X_agespace_mmc, X_spacetime_mmc, X_fixed_tmc, X_age_tmc, X_space_tmc,
    X_agespace_tmc,

    Q_space,

    u_fixed_mmc, 
    u_fixed_tmc, u_age_mmc,
    u_age_tmc, u_time_mmc, u_space_mmc, 
    u_space_tmc, 
    
    u_agetime_mmc, u_agespace_mmc,
    u_spacetime_mmc, u_agespace_tmc,
    
    logsigma_age_mmc, logsigma_time_mmc, logsigma_space_mmc,
    logsigma_agetime_mmc, logsigma_agespace_mmc, logsigma_spacetime_mmc,
    logsigma_age_tmc, logsigma_space_tmc, logsigma_agespace_tmc,
    
    logitrho_mmc_time1, logitrho_mmc_time2, logitrho_mmc_time3,
    logitrho_mmc_age1, logitrho_mmc_age2, logitrho_mmc_age3, logitrho_tmc_age1,
    logitrho_tmc_age2,
    
    report_vals
  );
  
  ///////////////////////////
  /// Reporting variables ///
  ///////////////////////////
  
  REPORT(report_vals.haz_mmc);     // Medical hazard rate
  REPORT(report_vals.haz_tmc);     // Traditional hazard rate
  REPORT(report_vals.haz);         // Total hazard rate
  REPORT(report_vals.inc_tmc);     // Traditional circumcision incidence rate
  REPORT(report_vals.inc_mmc);     // Medical circumcision incidence rate
  REPORT(report_vals.inc);         // Total circumcision incidence rate
  REPORT(report_vals.cum_inc_tmc); // Traditional circumcision cumulative incidence rate
  REPORT(report_vals.cum_inc_mmc); // Medical circumcision cumulative incidence rate
  REPORT(report_vals.cum_inc);     // Total circumcision cumulative incidence rate
  REPORT(report_vals.surv);        // Survival probabilities
    
  /////////////////////////////////////////
  /// Returning negative log likelihood ///
  /////////////////////////////////////////
  return(nll);
}
