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
  // private:
  public:
    // negative log likelihood
    Type nll; 
    // report values (hazard rates, incidence and cumulative incidence)
    vector<Type> haz_mmc;     // Medical hazard rate
    vector<Type> haz_tmc;     // Traditional hazard rate
    vector<Type> haz;         // Total hazard rate
    vector<Type> inc_tmc;     // Traditional circumcision incidence rate
    vector<Type> inc_mmc;     // Medical circumcision incidence rate
    vector<Type> inc;         // Total circumcision incidence rate
    vector<Type> cum_inc_tmc; // Traditional circumcision cumulative incidence rate
    vector<Type> cum_inc_mmc; // Medical circumcision cumulative incidence rate
    vector<Type> cum_inc;     // Total circumcision cumulative incidence rate
    vector<Type> surv;        // Survival probabilities
    vector<Type> leftcens;    // Left censoring probabilities (only used in likelihood)

    // also add report values here (??)

  // public:
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
    void rand_eff_space_p(SparseMatrix<Type> Q_space,
                          vector<Type> u_space_mmc,
                          vector<Type> u_space_tmc,
                          Type logsigma_space_mmc,
                          Type sigma_space_mmc,
                          Type logsigma_space_tmc,
                          Type sigma_space_tmc) {

      // Gaussian markov random field with prespecified precision matrix
      nll += GMRF(Q_space)(u_space_mmc);
      nll += GMRF(Q_space)(u_space_tmc);

      // Sum to zero constraints
      nll -= dnorm(u_space_mmc.sum(), Type(0), Type(0.001) * u_space_mmc.size(), TRUE);
      nll -= dnorm(u_space_tmc.sum(), Type(0), Type(0.001) * u_space_tmc.size(), TRUE);

      // Prior on the standard deviation for the spatial random effects
      nll -= dexp(sigma_space_mmc, Type(1), TRUE) + logsigma_space_mmc;
      nll -= dexp(sigma_space_tmc, Type(1), TRUE) + logsigma_space_tmc;
    };

    // Prior on the interaction random effects
    void rand_eff_interact_p(SparseMatrix<Type> Q_space,
                             array<Type> u_agespace_mmc,
                             array<Type> u_agespace_tmc,
                             array<Type> u_agetime_mmc,
                             array<Type> u_spacetime_mmc,
                             Type logsigma_agespace_mmc,
                             Type sigma_agespace_mmc,
                             Type logsigma_agespace_tmc,
                             Type sigma_agespace_tmc,
                             Type logsigma_agetime_mmc,
                             Type sigma_agetime_mmc,
                             Type logsigma_spacetime_mmc,
                             Type sigma_spacetime_mmc,
                             Type logitrho_mmc_age2,
                             Type rho_mmc_age2,
                             Type logitrho_tmc_age2,
                             Type rho_tmc_age2,
                             Type logitrho_mmc_age3,
                             Type rho_mmc_age3,
                             Type logitrho_mmc_time2,
                             Type rho_mmc_time2,
                             Type logitrho_mmc_time3,
                             Type rho_mmc_time3) {

      // Interactions: space-time (GMRF x AR1), age-time (AR1 x AR1) and age-space (AR1 x GMRF)
      nll += SEPARABLE(AR1(rho_mmc_time2), AR1(rho_mmc_age2))(u_agetime_mmc);
      nll += SEPARABLE(GMRF(Q_space), AR1(rho_mmc_age3))(u_agespace_mmc);
      nll += SEPARABLE(GMRF(Q_space), AR1(rho_mmc_time3))(u_spacetime_mmc);
      nll += SEPARABLE(GMRF(Q_space), AR1(rho_tmc_age2))(u_agespace_tmc);
      
      // Sum-to-zero constraints
      for (int i = 0; i < u_agespace_mmc.cols(); i++) {
        nll -= dnorm(u_agespace_mmc.col(i).sum(),
                     Type(0),
                     Type(0.001) * u_agespace_mmc.col(i).size(),
                     true);
      } 
      for (int i = 0; i < u_agetime_mmc.cols(); i++) {
        nll -= dnorm(u_agetime_mmc.col(i).sum(),
                     Type(0),
                     Type(0.001) * u_agetime_mmc.col(i).size(),
                     true);
      }  
      for (int i = 0; i < u_spacetime_mmc.cols(); i++) {
        nll -= dnorm(u_spacetime_mmc.col(i).sum(),
                     Type(0),
                     Type(0.001) * u_spacetime_mmc.col(i).size(),
                     true);
      }  
      for (int i = 0; i < u_agespace_tmc.cols(); i++) {
        nll -= dnorm(u_agespace_tmc.col(i).sum(),
                     Type(0),
                     Type(0.001) * u_agespace_tmc.col(i).size(),
                     true);
      }
      
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
    };

    // Function to calculate report values
    // TODO: This will change depending on whether type information is included
    void calc_report_vals(SparseMatrix<Type> X_fixed_mmc, 
                          SparseMatrix<Type> X_time_mmc,
                          SparseMatrix<Type> X_age_mmc, 
                          SparseMatrix<Type> X_space_mmc,
                          SparseMatrix<Type> X_agetime_mmc, 
                          SparseMatrix<Type> X_agespace_mmc,
                          SparseMatrix<Type> X_spacetime_mmc, 
                          // TMC terms
                          SparseMatrix<Type> X_fixed_tmc,
                          SparseMatrix<Type> X_age_tmc, 
                          SparseMatrix<Type> X_space_tmc,
                          SparseMatrix<Type> X_agespace_tmc,
                          // Integration matrices
                          density::SparseMatrix<Type> IntMat1,
                          density::SparseMatrix<Type> IntMat2,
                          // parameters
                          vector<Type> u_fixed_mmc, 
                          vector<Type> u_fixed_tmc,
                          vector<Type> u_age_mmc,
                          vector<Type> u_age_tmc,
                          vector<Type> u_time_mmc,
                          vector<Type> u_space_mmc, 
                          vector<Type> u_space_tmc,
                          array<Type> u_agetime_mmc,
                          array<Type> u_agespace_mmc,
                          array<Type> u_spacetime_mmc,
                          array<Type> u_agespace_tmc,
                          Type sigma_age_mmc,
                          Type sigma_time_mmc,
                          Type sigma_space_mmc,
                          Type sigma_agetime_mmc,
                          Type sigma_agespace_mmc,
                          Type sigma_spacetime_mmc,
                          Type sigma_age_tmc,
                          Type sigma_space_tmc,
                          Type sigma_agespace_tmc) {
       
      // Vector the interaction terms
      vector<Type> u_agespace_mmc_v(u_agespace_mmc);
      vector<Type> u_agetime_mmc_v(u_agetime_mmc);
      vector<Type> u_spacetime_mmc_v(u_spacetime_mmc);
      vector<Type> u_agespace_tmc_v(u_agespace_tmc);

      /// Estimate hazard rate ///
      /// TODO: Break this down into functions as well!
      // Medical hazard rate
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
      surv = exp(IntMat1 * logprob);
      vector<Type> surv_lag = exp(IntMat2 * logprob);
      leftcens = Type(1.0) - surv;

      // Incidence
      inc_tmc = haz_tmc * surv_lag;
      inc_mmc = haz_mmc * surv_lag;
      inc = haz * surv_lag;

      // Cumulative incidence
      cum_inc_tmc = IntMat1 * inc_tmc;
      cum_inc_mmc = IntMat1 * inc_mmc;
      cum_inc = cum_inc_tmc + cum_inc_mmc;
    };

    /// Calculate likelihood for types of circumcision ///
    void likelihood(SparseMatrix<Type> A_mmc,
                    SparseMatrix<Type> A_tmc,
                    SparseMatrix<Type> A_mc,
                    SparseMatrix<Type> B,
                    SparseMatrix<Type> C) {

      // likelihood for those medically circumcised
      nll -= (A_mmc * log(inc_mmc)).sum();

      // likelihood for those traditionally circumcised
      nll -= (A_tmc * log(inc_tmc)).sum();

      // // likelihood for those circumcised of unknown type
      nll -= (A_mc * log(inc)).sum();

      // likelihood for those right censored
      nll -= (B * log(surv)).sum();

      // // likelihood for those left censored
      nll -= (C * log(leftcens)).sum();
    };

    //// Getter Functions ////
   
    // getter for nll;
    Type get_nll() {
      return nll;
    };
    // getters for report vals
    vector<Type> get_haz_mmc() {
      return haz_mmc;
    };
    vector<Type> get_haz_tmc() {
      return haz_tmc;
    };
    vector<Type> get_haz() {
      return haz;
    };
    vector<Type> get_inc_mmc() {
      return inc_mmc;
    };
    vector<Type> get_inc_tmc() {
      return inc_tmc;
    };
    vector<Type> get_inc() {
      return inc;
    };
    vector<Type> get_cum_inc_mmc() {
      return cum_inc_mmc;
    };
    vector<Type> get_cum_inc_tmc() {
      return cum_inc_tmc;
    };
    vector<Type> get_cum_inc() {
      return cum_inc;
    };
    vector<Type> get_surv() {
      return surv;
    };
};

#endif

