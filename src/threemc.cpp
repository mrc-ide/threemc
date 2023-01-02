/// @file threemc.cpp


#define TMB_LIB_INIT R_init_threemc
#include <TMB.hpp>
#include "utils.h"
#include "implementation.cpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace density;

  // Negative log likelihood definition
  Type nll = Type(0);

  // Report values, to be passed by reference to function
  vector<Type> haz_mmc;
  vector<Type> haz_tmc;
  vector<Type> haz;
  vector<Type> inc_mmc;
  vector<Type> inc_tmc;
  vector<Type> inc;
  vector<Type> cum_inc_mmc;
  vector<Type> cum_inc_tmc;
  vector<Type> cum_inc;
  vector<Type> surv;

  // indicators
  DATA_INTEGER(is_type); // Does the data for the modelled country include type information
  DATA_INTEGER(paed_age_cutoff); // Fit with fixed paediatric MMC rate?

  if (is_type == 1) {

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

    /////////////////////////////////////
    // calculate nll and report values //
    /////////////////////////////////////
  
    // calculate nll when paed_age_cutoff is specified
    if (paed_age_cutoff == 1) {

      /////////////////////////////////
      /// Paed MMC Data definitions ///
      /////////////////////////////////
    
      // Design matrices 
      DATA_SPARSE_MATRIX(X_fixed_mmc_paed); // Design matrix for the fixed effects in the medical circumcision hazard rate (under specified age)
      DATA_SPARSE_MATRIX(X_age_mmc_paed); // Design matrix for the stratification random effects in "" "" ""
      DATA_SPARSE_MATRIX(X_space_mmc_paed); // Design matrix for the stratification random effects in "" "" ""
      DATA_SPARSE_MATRIX(X_agespace_mmc_paed); // Design matrix for the interaction random effects in "" "" ""
      
      ///////////////////////////
      /// Paed MMC Parameters ///
      ///////////////////////////

      // Fixed Effects
      PARAMETER_VECTOR(u_fixed_mmc_paed);

      // Age random effect
      PARAMETER_VECTOR(u_age_mmc_paed); 
    
      // Spatial random effects
      PARAMETER_VECTOR(u_space_mmc_paed);
    
      // Interactions
      PARAMETER_ARRAY(u_agespace_mmc_paed);
    
      // Standard deviations 
      PARAMETER(logsigma_age_mmc_paed)   // Type sigma_age_mmc_paed  = exp(logsigma_age_mmc_paed);
      PARAMETER(logsigma_space_mmc_paed); // Type sigma_space_mmc_paed = exp(logsigma_space_mmc_paed);
      PARAMETER(logsigma_agespace_mmc_paed); // Type sigma_agespace_mmc_paed = exp(logsigma_agespace_mmc_paed);

      // Autocorrelation parameters 
      PARAMETER(logitrho_mmc_paed_age1);   // Type rho_mmc_age1   = geninvlogit(logitrho_mmc_age1,  Type(-1.0), Type(1.0));
      PARAMETER(logitrho_mmc_paed_age2);   // Type rho_mmc_age2   = geninvlogit(logitrho_mmc_age2,  Type(-1.0), Type(1.0));
      
      nll = threemc(
        
        A_mmc, A_tmc, A_mc, B, C, IntMat1, IntMat2, 
        
        X_fixed_mmc, X_time_mmc, X_age_mmc, X_space_mmc, X_agetime_mmc,
        X_agespace_mmc, X_spacetime_mmc,

        X_fixed_mmc_paed, X_age_mmc_paed, X_space_mmc_paed, X_agespace_mmc_paed,
        
        X_fixed_tmc, X_age_tmc, X_space_tmc,
        X_agespace_tmc,

        Q_space,

        u_fixed_mmc, u_fixed_mmc_paed, u_fixed_tmc,
        u_age_mmc, u_age_mmc_paed, u_age_tmc,
        u_time_mmc,
        u_space_mmc, u_space_mmc_paed,
        u_space_tmc, 
        
        u_agetime_mmc,
        u_agespace_mmc, u_agespace_mmc_paed,
        u_spacetime_mmc, u_agespace_tmc,

        logsigma_age_mmc, logsigma_time_mmc, logsigma_space_mmc,
        logsigma_age_mmc_paed, logsigma_space_mmc_paed,
        logsigma_agetime_mmc, logsigma_agespace_mmc, logsigma_spacetime_mmc,
        logsigma_agespace_mmc_paed,
        logsigma_age_tmc, logsigma_space_tmc, logsigma_agespace_tmc,
        
        logitrho_mmc_time1, logitrho_mmc_time2, logitrho_mmc_time3,
        logitrho_mmc_age1,  logitrho_mmc_paed_age1,
        logitrho_mmc_age2,  logitrho_mmc_paed_age2,
        logitrho_mmc_age3,
        logitrho_tmc_age1, logitrho_tmc_age2,
        
        // report vals
        haz_mmc, haz_tmc, haz, inc_mmc, inc_tmc, inc,
        cum_inc_mmc, cum_inc_tmc, cum_inc, surv
       );
    } else {
      // calculate nll when no paed_age_cutoff is specified
      nll = threemc(

         A_mmc, A_tmc, A_mc, B, C, IntMat1, IntMat2, 
         
         X_fixed_mmc, X_time_mmc, X_age_mmc, X_space_mmc, X_agetime_mmc,
         X_agespace_mmc, X_spacetime_mmc,

         X_fixed_tmc, X_age_tmc, X_space_tmc,
         X_agespace_tmc,

         Q_space,

         u_fixed_mmc, u_fixed_tmc,
         u_age_mmc, u_age_tmc,
         u_time_mmc,
         u_space_mmc, 
         u_space_tmc, 
         
         u_agetime_mmc,
         u_agespace_mmc, 
         u_spacetime_mmc, u_agespace_tmc,

         logsigma_age_mmc, logsigma_time_mmc, logsigma_space_mmc,
         logsigma_agetime_mmc, logsigma_agespace_mmc, logsigma_spacetime_mmc,
         logsigma_age_tmc, logsigma_space_tmc, logsigma_agespace_tmc,
         
         logitrho_mmc_time1, logitrho_mmc_time2, logitrho_mmc_time3,
         logitrho_mmc_age1,
         logitrho_mmc_age2,
         logitrho_mmc_age3,
         logitrho_tmc_age1, logitrho_tmc_age2,
         
         // report vals
         haz_mmc, haz_tmc, haz, inc_mmc, inc_tmc, inc,
         cum_inc_mmc, cum_inc_tmc, cum_inc, surv
      );
    }
    
    ///////////////////////////
    /// Reporting variables ///
    ///////////////////////////

    REPORT(haz_mmc);     // Medical hazard rate
    REPORT(haz_tmc);     // Traditional hazard rate
    REPORT(haz);         // Total hazard rate
    REPORT(inc_tmc);     // Traditional circumcision incidence rate
    REPORT(inc_mmc);     // Medical circumcision incidence rate
    REPORT(inc);         // Total circumcision incidence rate
    REPORT(cum_inc_tmc); // Traditional circumcision cumulative incidence rate
    REPORT(cum_inc_mmc); // Medical circumcision cumulative incidence rate
    REPORT(cum_inc);     // Total circumcision cumulative incidence rate
    REPORT(surv);        // Survival probabilities

  } else if (is_type == 0) {

    ////////////////////////
    /// Data definitions ///
    ////////////////////////
    // Survival analysis matrices
    DATA_SPARSE_MATRIX(A); // Matrix selecting instantaneous hazard for unknown circumcised pop
    DATA_SPARSE_MATRIX(B); // Matrix selecting relevant cumulative hazard entry for observed and right censored pop
    DATA_SPARSE_MATRIX(C); // Matrix selecting relevant cumulative hazard entry for interval censored pop
    DATA_SPARSE_MATRIX(IntMat1); // Integration matrix for cumulative hazard 
    DATA_SPARSE_MATRIX(IntMat2); // Integration matrix for lagged cumulative hazard 
    
    // Design matrices 
    DATA_SPARSE_MATRIX(X_fixed); // Design matrix for the fixed effects in the medical circumcision hazard rate
    DATA_SPARSE_MATRIX(X_time); // Design matrix for the temporal random effects in the medical circumcision hazard rate
    DATA_SPARSE_MATRIX(X_age); // Design matrix for the stratification random effects in the medical circumcision hazard rate
    DATA_SPARSE_MATRIX(X_space); // Design matrix for the stratification random effects in the medical circumcision hazard rate
    DATA_SPARSE_MATRIX(X_agetime); // Design matrix for the interaction random effects in the medical circumcision hazard rate
    DATA_SPARSE_MATRIX(X_agespace); // Design matrix for the interaction random effects in the medical circumcision hazard rate
    DATA_SPARSE_MATRIX(X_spacetime); // Design matrix for the interaction random effects in the medical circumcision hazard rate
    
    // Precision matrices 
    DATA_SPARSE_MATRIX(Q_space); // Aggregation matrix for number of circumcisions performed
    
    //////////////////
    /// Parameters ///
    //////////////////
    // Fixed Effects
    PARAMETER_VECTOR(u_fixed);
    
    // Age random effect
    PARAMETER_VECTOR(u_age); 
    
    // Temporal random effects 
    PARAMETER_VECTOR(u_time);
    
    // Spatial random effects
    PARAMETER_VECTOR(u_space);
    
    // Interactions
    PARAMETER_ARRAY(u_agetime);
    PARAMETER_ARRAY(u_agespace);
    PARAMETER_ARRAY(u_spacetime);
    
    // Standard deviations 
    PARAMETER(logsigma_age);
    PARAMETER(logsigma_time);
    PARAMETER(logsigma_space);
    PARAMETER(logsigma_agetime);
    PARAMETER(logsigma_agespace);
    PARAMETER(logsigma_spacetime);
    
    // Autocorrelation parameters 
    PARAMETER(logitrho_time1);
    PARAMETER(logitrho_time2);
    PARAMETER(logitrho_time3);
    PARAMETER(logitrho_age1);
    PARAMETER(logitrho_age2);
    PARAMETER(logitrho_age3);
 
    // Negative log likelihood definition
    // Type nll = Type(0);
    
    // Report values, to be passed by reference to function
    // vector<Type> haz;
    // vector<Type> inc;
    // vector<Type> cum_inc;
    // vector<Type> surv;


    nll = threemc(

      A, B, C, IntMat1, IntMat2, 
      
      X_fixed, X_time, X_age, X_space, X_agetime,
      X_agespace, X_spacetime, 

      Q_space,

      u_fixed, u_age, u_time, u_space, 
      
      u_agetime, u_agespace, u_spacetime,
      
      logsigma_age, logsigma_time, logsigma_space,
      logsigma_agetime, logsigma_agespace, logsigma_spacetime,
      
      logitrho_time1, logitrho_time2, logitrho_time3,
      logitrho_age1, logitrho_age2, logitrho_age3,
      
      // report vals
      haz, inc, cum_inc, surv
    );

    ///////////////////////////
    /// Reporting variables ///
    ///////////////////////////
    REPORT(haz);         // Total hazard rate
    REPORT(inc);         // Total circumcision incidence rate
    REPORT(cum_inc);     // Total circumcision cumulative incidence rate
    REPORT(surv);        // Survival probabilities
  }
 
  /////////////////////////////////////////
  /// Returning negative log likelihood ///
  /////////////////////////////////////////
  return(nll);
}
