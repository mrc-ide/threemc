#include <TMB.hpp>
// #include "utils.h"

/* Utility Functions (perhaps move to the end? Or their own header?) */



// /// @file utils.h
// #pragma once

/***************************************************/
/* Function to inverse logit transform for vectors */
/***************************************************/
template <class vector>
vector invlogit_vec(vector x){
  vector y = 1.0 / (1.0 + exp(-x));
  return y;
}

/*******************************************************************/
/* Function to inverse logit transform on a general interval [a,b] */
/*******************************************************************/
template <class Type>
Type geninvlogit(Type x, Type a, Type b){
  Type y; 
  y = 1.0 / (1.0 + exp(-x));
  y = y * (b - a) + a;
  return y;
}

/*******************************************************************/
/*   Class to calculate negative log likelihood in threemc model   */
/*******************************************************************/

#ifndef THREEMCHEADERDEF
#define THREEMCHEADERDEF

using namespace density;

// TODO: Move implementation to separate file
template <class Type>
class Threemc {
  private:
    Type nll; // negative log likelihood
    // also add report values here (??)

  public:
    // Default Constructor
    Threemc() {
      Type nll = Type(0); // initialise nll to 0
      // TODO: also intitialise report values here
    };

    // TODO: Write constructor with more arguments?

    // Default Destructor (needed??)
    // ~Threemc() {
    // };

    // Prior on fixed effects (won't need second version until implementing no type version)
    void fix_eff_p(vector<Type> u_fixed_mmc, vector<Type> u_fixed_tmc) {
      // fixed effects for the medical circumcision rate
      nll -= dnorm(u_fixed_mmc, Type(0), Type(5), true).sum();
      // fixed effects for the traditional circumcision rate
      nll -= dnorm(u_fixed_tmc, Type(0), Type(5), true).sum();
    };

    // Prior on temporal random effects
    // TODO: This will become a wrapper (?) for other functions when other models are implemented
    void rand_eff_time_p(vector<Type> u_time_mmc,
                         Type logsigma_time_mmc,
                         Type sigma_time_mmc,
                         Type logitrho_mmc_time1,
                         Type rho_mmc_time1) {

      // AR1 Process
      nll += AR1(rho_mmc_time1)(u_time_mmc);

      // Sum to zero constraint
      nll -= dnorm(u_time_mmc.sum(), Type(0), Type(0.001) * u_time_mmc.size(), TRUE);

      // Prior on the standard deviation for the temporal random effects
      nll -= dexp(sigma_time_mmc, Type(1), TRUE) + logsigma_time_mmc;

      // Prior on the logit autocorrelation parameters
      nll -= dnorm(logitrho_mmc_time1, Type(3), Type(3), TRUE);
    };


    // Prior on the age random effects
    // TODO: Will change in model where there is a paediatric-adult age split
    void rand_eff_age_p(vector<Type> u_age_mmc,
                        vector<Type> u_age_tmc,
                        Type logsigma_age_mmc,
                        Type sigma_age_mmc,
                        Type logsigma_age_tmc,
                        Type sigma_age_tmc,
                        Type logitrho_mmc_age1,
                        Type rho_mmc_age1,
                        Type logitrho_tmc_age1,
                        Type rho_tmc_age1) {

      // AR1 processes
      nll += AR1(rho_mmc_age1)(u_age_mmc);
      nll += AR1(rho_tmc_age1)(u_age_tmc);

      // sum to zero constraint
      nll -= dnorm(u_age_mmc.sum(), Type(0), Type(0.001) * u_age_mmc.size(), true);
      nll -= dnorm(u_age_tmc.sum(), Type(0), Type(0.001) * u_age_tmc.size(), true);

      // prior on the standard deviation for the age random effects
      nll -= dexp(sigma_age_mmc, Type(1), true) + logsigma_age_mmc;
      nll -= dexp(sigma_age_tmc, Type(1), true) + logsigma_age_tmc;

      // prior on the logit autocorrelation parameters
      nll -= dnorm(logitrho_mmc_age1, Type(3), Type(2), true);
      nll -= dnorm(logitrho_tmc_age1, Type(3), Type(2), true);
    };

    // Prior on the spatial random effects
    // void rand_eff_space_p(SparseMatrix<Type> Q_space,
    //                       vector<Type> u_space_mmc,
    //                       vector<Type> u_space_tmc,
    //                       Type logsigma_space_mmc,
    //                       Type sigma_space_mmc,
    //                       Type logsigma_space_tmc,
    //                       Type sigma_space_tmc) {

    //   // Gaussian markov random field with prespecified precision matrix
    //   nll += GMRF(Q_space)(u_space_mmc);
    //   nll += GMRF(Q_space)(u_space_tmc);

    //   // Sum to zero constraints
    //   nll -= dnorm(u_space_mmc.sum(), Type(0), Type(0.001) * u_space_mmc.size(), TRUE);
    //   nll -= dnorm(u_space_tmc.sum(), Type(0), Type(0.001) * u_space_tmc.size(), TRUE);

    //   // Prior on the standard deviation for the spatial random effects
    //   nll -= dexp(sigma_space_mmc, Type(1), TRUE) + logsigma_space_mmc;
    //   nll -= dexp(sigma_space_tmc, Type(1), TRUE) + logsigma_space_tmc;
    // };

    // // Prior on the interaction random effects
    // void rand_eff_interact_p(SparseMatrix<Type> Q_space,
    //                          array<Type> u_agespace_mmc,
    //                          array<Type> u_agespace_tmc,
    //                          array<Type> u_agetime_mmc,
    //                          array<Type> u_spacetime_mmc,
    //                          Type logsigma_agespace_mmc,
    //                          Type sigma_agespace_mmc,
    //                          Type logsigma_agespace_tmc,
    //                          Type sigma_agespace_tmc,
    //                          Type logsigma_agetime_mmc,
    //                          Type sigma_agetime_mmc,
    //                          Type logsigma_spacetime_mmc,
    //                          Type sigma_spacetime_mmc,
    //                          Type logitrho_mmc_age2,
    //                          Type rho_mmc_age2,
    //                          Type logitrho_tmc_age2,
    //                          Type rho_tmc_age2,
    //                          Type logitrho_mmc_age3,
    //                          Type rho_mmc_age3,
    //                          Type logitrho_mmc_time2,
    //                          Type rho_mmc_time2,
    //                          Type logitrho_mmc_time3,
    //                          Type rho_mmc_time3) {

    //   // Interactions: space-time (GMRF x AR1), age-time (AR1 x AR1) and age-space (AR1 x GMRF)
    //   nll += SEPARABLE(AR1(rho_mmc_time2), AR1(rho_mmc_age2))(u_agetime_mmc);
    //   nll += SEPARABLE(GMRF(Q_space), AR1(rho_mmc_age3))(u_agespace_mmc);
    //   nll += SEPARABLE(GMRF(Q_space), AR1(rho_mmc_time3))(u_spacetime_mmc);
    //   nll += SEPARABLE(GMRF(Q_space), AR1(rho_tmc_age2))(u_agespace_tmc);
    //   
    //   // Sum-to-zero constraints
    //   for (int i = 0; i < u_agespace_mmc.cols(); i++) {
    //     nll -= dnorm(u_agespace_mmc.col(i).sum(),
    //                  Type(0),
    //                  Type(0.001) * u_agespace_mmc.col(i).size(),
    //                  true);
    //   } 
    //   for (int i = 0; i < u_agetime_mmc.cols(); i++) {
    //     nll -= dnorm(u_agetime_mmc.col(i).sum(),
    //                  Type(0),
    //                  Type(0.001) * u_agetime_mmc.col(i).size(),
    //                  true);
    //   }  
    //   for (int i = 0; i < u_spacetime_mmc.cols(); i++) {
    //     nll -= dnorm(u_spacetime_mmc.col(i).sum(),
    //                  Type(0),
    //                  Type(0.001) * u_spacetime_mmc.col(i).size(),
    //                  true);
    //   }  
    //   for (int i = 0; i < u_agespace_tmc.cols(); i++) {
    //     nll -= dnorm(u_agespace_tmc.col(i).sum(),
    //                  Type(0),
    //                  Type(0.001) * u_agespace_tmc.col(i).size(),
    //                  true);
    //   }
    //   
    //   // Vectorising the interaction
    //   // TODO: Move this to function where report values are calculated
    //   // vector<Type> u_agespace_mmc_v(u_agespace_mmc);
    //   // vector<Type> u_agetime_mmc_v(u_agetime_mmc);
    //   // vector<Type> u_spacetime_mmc_v(u_spacetime_mmc);
    //   // vector<Type> u_agespace_tmc_v(u_agespace_tmc);

    //   // Prior on the standard deviation for the interaction random effects
    //   nll -= dexp(sigma_agespace_mmc,  Type(1), TRUE) + logsigma_agespace_mmc;
    //   nll -= dexp(sigma_agetime_mmc,   Type(1), TRUE) + logsigma_agetime_mmc;
    //   nll -= dexp(sigma_spacetime_mmc, Type(1), TRUE) + logsigma_spacetime_mmc;
    //   nll -= dexp(sigma_agespace_tmc,  Type(1), TRUE) + logsigma_agespace_tmc;

    //   // Prior on the logit autocorrelation parameters
    //   nll -= dnorm(logitrho_mmc_time2, Type(3), Type(2), TRUE);
    //   nll -= dnorm(logitrho_mmc_age2,  Type(3), Type(2), TRUE);
    //   nll -= dnorm(logitrho_mmc_time3, Type(3), Type(2), TRUE);
    //   nll -= dnorm(logitrho_mmc_age3,  Type(3), Type(2), TRUE);
    //   nll -= dnorm(logitrho_tmc_age2,  Type(3), Type(2), TRUE);
    // };

    // Function to calculate report values
    // TODO: This will change depending on whether type information is included
    // calc_report_vals() {
    //   
    // };

    // Function
    // likelihood() {};

    // getter for nll;
    Type get_nll() {
      return nll;
    };

    // Report values (
    // TOOD: Need two versions depending on whether MC is split by type)
    // void report() {
    //   REPORT(haz_mmc);     // Medical hazard rate
    //   REPORT(haz_tmc);     // Traditional hazard rate
    //   REPORT(haz);         // Total hazard rate
    //   REPORT(inc_tmc);     // Traditional circumcision incidence rate
    //   REPORT(inc_mmc);     // Medical circumcision incidence rate
    //   REPORT(inc);         // Total circumcision incidence rate
    //   REPORT(cum_inc_tmc); // Traditional circumcision cumulative incidence rate
    //   REPORT(cum_inc_mmc); // Medical circumcision cumulative incidence rate
    //   REPORT(cum_inc);     // Total circumcision cumulative incidence rate
    //   REPORT(surv);        // Survival probabilities
    // };
};

#endif

