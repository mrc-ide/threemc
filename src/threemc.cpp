#define TMB_LIB_INIT R_init_threemc
#include "threemc.h"
// #include "implementation.cpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace density;

  //////////////////////////
  //// Data definitions ////
  //////////////////////////

  // Define struct containing data matrices
  DATA_STRUCT(threemc_data, Threemc_data);
  
  ////////////////////
  //// Parameters ////
  ////////////////////

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
  PARAMETER(logsigma_age_mmc);       Type sigma_age_mmc       = exp(logsigma_age_mmc);
  PARAMETER(logsigma_time_mmc);      Type sigma_time_mmc      = exp(logsigma_time_mmc);
  PARAMETER(logsigma_space_mmc);     Type sigma_space_mmc     = exp(logsigma_space_mmc);
  PARAMETER(logsigma_agetime_mmc);   Type sigma_agetime_mmc   = exp(logsigma_agetime_mmc);
  PARAMETER(logsigma_agespace_mmc);  Type sigma_agespace_mmc  = exp(logsigma_agespace_mmc);
  PARAMETER(logsigma_spacetime_mmc); Type sigma_spacetime_mmc = exp(logsigma_spacetime_mmc);
  PARAMETER(logsigma_age_tmc);       Type sigma_age_tmc       = exp(logsigma_age_tmc);
  PARAMETER(logsigma_space_tmc);     Type sigma_space_tmc     = exp(logsigma_space_tmc);
  PARAMETER(logsigma_agespace_tmc);  Type sigma_agespace_tmc  = exp(logsigma_agespace_tmc);

  // Autocorrelation parameters 
  PARAMETER(logitrho_mmc_time1);
  Type rho_mmc_time1  = geninvlogit(logitrho_mmc_time1, Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_time2);
  Type rho_mmc_time2  = geninvlogit(logitrho_mmc_time2, Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_time3);
  Type rho_mmc_time3  = geninvlogit(logitrho_mmc_time3, Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_age1);
  Type rho_mmc_age1   = geninvlogit(logitrho_mmc_age1,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_age2);
  Type rho_mmc_age2   = geninvlogit(logitrho_mmc_age2,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_mmc_age3);
  Type rho_mmc_age3   = geninvlogit(logitrho_mmc_age3,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_tmc_age1);
  Type rho_tmc_age1   = geninvlogit(logitrho_tmc_age1,  Type(-1.0), Type(1.0));
  PARAMETER(logitrho_tmc_age2);
  Type rho_tmc_age2   = geninvlogit(logitrho_tmc_age2,  Type(-1.0), Type(1.0));


  //// Calculate nll and report values ////

  // define object of class Threemc, which will store our negative log likelihood
  Threemc<Type> threemc;

  //// Priors ////
  
  // Apply prior on fixed effects
  threemc.fix_eff_p(u_fixed_mmc, u_fixed_tmc);

  // prior on temporal random effects
  threemc.rand_eff_time_p(u_time_mmc,
                          logsigma_time_mmc,
                          sigma_time_mmc,
                          logitrho_mmc_time1,
                          rho_mmc_time1);

  // prior on the age random effects
  threemc.rand_eff_age_p(u_age_mmc,
                         u_age_tmc,
                         logsigma_age_mmc,
                         sigma_age_mmc,
                         logsigma_age_tmc,
                         sigma_age_tmc,
                         logitrho_mmc_age1,
                         rho_mmc_age1,
                         logitrho_tmc_age1,
                         rho_tmc_age1);

  // prior on spatial random effects
  threemc.rand_eff_space_p(threemc_data.Q_space,
                           u_space_mmc,
                           u_space_tmc,
                           logsigma_space_mmc,
                           sigma_space_mmc,
                           logsigma_space_tmc,
                           sigma_space_tmc);

  // prior on interaction random effects
  threemc.rand_eff_interact_p(threemc_data.Q_space,
                              u_agespace_mmc,
                              u_agespace_tmc,
                              u_agetime_mmc,
                              u_spacetime_mmc,
                              logsigma_agespace_mmc,
                              sigma_agespace_mmc,
                              logsigma_agespace_tmc,
                              sigma_agespace_tmc,
                              logsigma_agetime_mmc,
                              sigma_agetime_mmc,
                              logsigma_spacetime_mmc,
                              sigma_spacetime_mmc,
                              logitrho_mmc_age2,
                              rho_mmc_age2,
                              logitrho_tmc_age2,
                              rho_tmc_age2,
                              logitrho_mmc_age3,
                              rho_mmc_age3,
                              logitrho_mmc_time2,
                              rho_mmc_time2,
                              logitrho_mmc_time3,
                              rho_mmc_time3);

  //// Calculate report values (hazard, (cumulative) incidence) ////

  // Calculate hazards
  threemc.calc_haz(threemc_data.X_fixed_mmc, 
                   threemc_data.X_time_mmc,
                   threemc_data.X_age_mmc, 
                   threemc_data.X_space_mmc,
                   threemc_data.X_agetime_mmc, 
                   threemc_data.X_agespace_mmc,
                   threemc_data.X_spacetime_mmc, 
                   u_fixed_mmc, 
                   u_age_mmc,
                   u_time_mmc,
                   u_space_mmc, 
                   u_agetime_mmc,
                   u_agespace_mmc,
                   u_spacetime_mmc,
                   sigma_age_mmc,
                   sigma_time_mmc,
                   sigma_space_mmc,
                   sigma_agetime_mmc,
                   sigma_agespace_mmc,
                   sigma_spacetime_mmc);


  threemc.calc_haz(threemc_data.X_fixed_tmc, 
                   threemc_data.X_age_tmc, 
                   threemc_data.X_space_tmc,
                   threemc_data.X_agespace_tmc,
                   u_fixed_tmc, 
                   u_age_tmc,
                   u_space_tmc, 
                   u_agespace_tmc,
                   sigma_age_tmc,
                   sigma_space_tmc,
                   sigma_agespace_tmc);

  threemc.calc_haz();

  // calculate survival probabilities
  threemc.calc_surv(threemc_data.IntMat1, threemc_data.IntMat2);

  // calculate incidences and cumulative incidences
  threemc.calc_inc(threemc_data.IntMat1, 1); // run where we have type

  //// Calculate likelihood ////
  threemc.likelihood(threemc_data.A_mmc, threemc.get_inc_mmc()); // medical circumcisions
  threemc.likelihood(threemc_data.A_tmc, threemc.get_inc_tmc()); // traditional circumcisions
  threemc.likelihood(threemc_data.A_mc, threemc.get_inc());      // circs of unknown type
  threemc.likelihood(threemc_data.B, threemc.get_surv());        // right censored (i.e. uncircumcised)
  threemc.likelihood(threemc_data.C, threemc.get_leftcens());    // left censored (i.e. unknown circ age)

  //// report hazard rates, incidence and cumulative incidence ////
  REPORT(threemc.get_haz_mmc());     // Medical hazard rate
  REPORT(threemc.get_haz_tmc());     // Traditional hazard rate
  REPORT(threemc.get_haz());         // Total hazard rate
  REPORT(threemc.get_inc_mmc());     // Medical circumcision incidence rate
  REPORT(threemc.get_inc_tmc());     // Traditional circumcision incidence rate
  REPORT(threemc.get_inc());         // Total circumcision incidence rate
  REPORT(threemc.get_cum_inc_mmc()); // Medical circumcision cumulative incidence rate
  REPORT(threemc.get_cum_inc_tmc()); // Traditional circumcision cumulative incidence rate
  REPORT(threemc.get_cum_inc());     // Total circumcision cumulative incidence rate
  REPORT(threemc.get_surv());        // Survival probabilities
  
  //// return nll ////
  return threemc.get_nll();
}
