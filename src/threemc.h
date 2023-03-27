#include <TMB.hpp>
// #include "utils.h"

/* Utility Functions (perhaps move to the end? Or their own header?) */

// TODO: Make wrapper functions for overloaded functions

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
    vector<Type> logprob;     // ?
    vector<Type> surv;        // Survival probabilities
    vector<Type> surv_lag;    // Lagged survival probabilities
    vector<Type> leftcens;    // left censored incidence rate

  public:
    // Default Constructor
    Threemc() {
      Type nll = Type(0); // initialise nll to 0
    };

    // TODO: Write constructor with more arguments, initialising report values (needed?)

    // Default Destructor (needed??)
    // ~Threemc() {
    // };

    // Prior on fixed effects 
    // IDEA: Could just do this twice for tmc and mmc if we have type, and once if we no type?!!
    void fix_eff_p(vector<Type> u_fixed) {
	    // Fixed effects for circumcision rate
      nll -= dnorm(u_fixed, Type(0), Type(5), true).sum();
    }

    // Prior on temporal random effects
    // This function works for both type and no type specifications (just with different inputs)
    void rand_eff_time_p(vector<Type> u_time,
                         Type logsigma_time,
                         Type sigma_time,
                         Type logitrho_time1,
                         Type rho_time1) {

      // AR1 Process
      nll += AR1(rho_time1)(u_time);

      // Sum to zero constraint
      nll -= dnorm(u_time.sum(), Type(0), Type(0.001) * u_time.size(), true);

      // Prior on the standard deviation for the temporal random effects
      nll -= dexp(sigma_time, Type(1), true) + logsigma_time;

      // Prior on the logit autocorrelation parameters
      nll -= dnorm(logitrho_time1, Type(3), Type(3), true);
    };


    // Prior on the age random effects
    // TODO: Will change in model where there is a paediatric-adult age split
    void rand_eff_age_p(vector<Type> u_age,
                        Type logsigma_age,
                        Type sigma_age,
                        Type logitrho_age1,
                        Type rho_age1) {

      // AR1 processes
      nll += AR1(rho_age1)(u_age);

      // sum to zero constraint
      nll -= dnorm(u_age.sum(), Type(0), Type(0.001) * u_age.size(), true);

      // prior on the standard deviation for the age random effects
      nll -= dexp(sigma_age, Type(1), true) + logsigma_age;

      // prior on the logit autocorrelation parameters
      nll -= dnorm(logitrho_age1, Type(3), Type(2), true);
    };

    // Prior on the spatial random effects
    void rand_eff_space_p(SparseMatrix<Type> Q_space,
                          vector<Type> u_space,
                          Type logsigma_space,
                          Type sigma_space) {

      // Gaussian markov random field with prespecified precision matrix
      nll += GMRF(Q_space)(u_space);

      // Sum to zero constraints
      nll -= dnorm(u_space.sum(), Type(0), Type(0.001) * u_space.size(), true);

      // Prior on the standard deviation for the spatial random effects
      nll -= dexp(sigma_space, Type(1), true) + logsigma_space;
    };

    // Prior on the interaction random effects for either MMC or MC (for model w/ no type)
    void rand_eff_interact_p(SparseMatrix<Type> Q_space,
                             array<Type> u_agespace,
                             array<Type> u_agetime,
                             array<Type> u_spacetime,
                             Type logsigma_agespace,
                             Type sigma_agespace,
                             Type logsigma_agetime,
                             Type sigma_agetime,
                             Type logsigma_spacetime,
                             Type sigma_spacetime,
                             Type logitrho_age2,
                             Type rho_age2,
                             Type logitrho_age3,
                             Type rho_age3,
                             Type logitrho_time2,
                             Type rho_time2,
                             Type logitrho_time3,
                             Type rho_time3) {

      // Interactions: space-time (GMRF x AR1), age-time (AR1 x AR1) and age-space (AR1 x GMRF)
      nll += SEPARABLE(AR1(rho_time2), AR1(rho_age2))(u_agetime);
      nll += SEPARABLE(GMRF(Q_space), AR1(rho_age3))(u_agespace);
      nll += SEPARABLE(GMRF(Q_space), AR1(rho_time3))(u_spacetime);
      
      // Sum-to-zero constraint
      // TODO: Iterate over map of these
      for (int i = 0; i < u_agespace.cols(); i++) {
        nll -= dnorm(u_agespace.col(i).sum(),
                     Type(0),
                     Type(0.001) * u_agespace.col(i).size(),
                     true);
      } 
      for (int i = 0; i < u_agetime.cols(); i++) {
        nll -= dnorm(u_agetime.col(i).sum(),
                     Type(0),
                     Type(0.001) * u_agetime.col(i).size(),
                     true);
      }  
      for (int i = 0; i < u_spacetime.cols(); i++) {
        nll -= dnorm(u_spacetime.col(i).sum(),
                     Type(0),
                     Type(0.001) * u_spacetime.col(i).size(),
                     true);
      }  
      
      // Prior on the standard deviation for the interaction random effects
      // TODO: Can rewrite this using a pointer map or something to iterate over
      nll -= dexp(sigma_agespace,  Type(1), true) + logsigma_agespace;
      nll -= dexp(sigma_agetime,   Type(1), true) + logsigma_agetime;
      nll -= dexp(sigma_spacetime, Type(1), true) + logsigma_spacetime;

      // Prior on the logit autocorrelation parameters
      // TODO: Can also iterate over these so we don't need to overload this function
      nll -= dnorm(logitrho_time2, Type(3), Type(2), true);
      nll -= dnorm(logitrho_age2,  Type(3), Type(2), true);
      nll -= dnorm(logitrho_time3, Type(3), Type(2), true);
      nll -= dnorm(logitrho_age3,  Type(3), Type(2), true);
    };

    // Overload prior on the interaction random effects for just TMC
    void rand_eff_interact_p(SparseMatrix<Type> Q_space,
                             array<Type> u_agespace_tmc,
                             Type logsigma_agespace_tmc,
                             Type sigma_agespace_tmc,
                             Type logitrho_tmc_age2,
                             Type rho_tmc_age2) {

      // Interactions: space-time (GMRF x AR1), age-time (AR1 x AR1) and age-space (AR1 x GMRF)
      nll += SEPARABLE(GMRF(Q_space), AR1(rho_tmc_age2))(u_agespace_tmc);
      
      // Sum-to-zero constraints
      for (int i = 0; i < u_agespace_tmc.cols(); i++) {
        nll -= dnorm(u_agespace_tmc.col(i).sum(),
                     Type(0),
                     Type(0.001) * u_agespace_tmc.col(i).size(),
                     true);
      }
      
      // Prior on the standard deviation for the interaction random effects
      nll -= dexp(sigma_agespace_tmc,  Type(1), true) + logsigma_agespace_tmc;

      // Prior on the logit autocorrelation parameters
      nll -= dnorm(logitrho_tmc_age2,  Type(3), Type(2), true);
    };


    // Function to calculate report values 
    // TODO: This will change depending on whether type information is included
    // Need to just overload this function
    // TODO: Can definitely design this function better to avoid repitition
    // For MMC: 
    void calc_haz(SparseMatrix<Type> X_fixed, 
                  SparseMatrix<Type> X_time,
                  SparseMatrix<Type> X_age, 
                  SparseMatrix<Type> X_space,
                  SparseMatrix<Type> X_agetime, 
                  SparseMatrix<Type> X_agespace,
                  SparseMatrix<Type> X_spacetime, 
                  // parameters
                  vector<Type> u_fixed, 
                  vector<Type> u_age,
                  vector<Type> u_time,
                  vector<Type> u_space, 
                  array<Type> u_agetime,
                  array<Type> u_agespace,
                  array<Type> u_spacetime,
                  Type sigma_age,
                  Type sigma_time,
                  Type sigma_space,
                  Type sigma_agetime,
                  Type sigma_agespace,
                  Type sigma_spacetime) {

      // Vector the interaction terms
      vector<Type> u_agespace_v(u_agespace);
      vector<Type> u_agetime_v(u_agetime);
      vector<Type> u_spacetime_v(u_spacetime);

      /// Estimate hazard rate ///
      /// TODO: Break this down into functions as well!
      // Medical hazard rate
      // TODO: How to change haz_mmc to haz, but effect haz_mmc member??
      haz_mmc = X_fixed * u_fixed +
        X_time * u_time * sigma_time +
        X_space * u_space * sigma_space +
        X_age * u_age * sigma_age +
        X_agetime * u_agetime_v * sigma_agetime +
        X_agespace * u_agespace_v * sigma_agespace +
        X_spacetime * u_spacetime_v * sigma_spacetime;

      // Rates on [0,1] scale
      haz_mmc = invlogit_vec(haz_mmc);

      // // Incidence
      // inc_mmc = haz_mmc * surv_lag;

      // // Cumulative incidence
      // cum_inc_mmc = IntMat1 * inc;
    };

    // For TMC: 
    void calc_haz(SparseMatrix<Type> X_fixed,
                  SparseMatrix<Type> X_age, 
                  SparseMatrix<Type> X_space,
                  SparseMatrix<Type> X_agespace,
                  // parameters
                  vector<Type> u_fixed,
                  vector<Type> u_age,
                  vector<Type> u_space,
                  array<Type> u_agespace,
                  Type sigma_age,
                  Type sigma_space,
                  Type sigma_agespace) {

      // Vectorise interaction terms
      vector<Type> u_agespace_v(u_agespace);

      /// Estimate hazard rate ///
      haz_tmc = X_fixed * u_fixed +
        X_space * u_space * sigma_space +
        X_age * u_age * sigma_age +
        X_agespace * u_agespace_v * sigma_agespace;

      // Rates on [0,1] scale
      haz_tmc = invlogit_vec(haz_tmc);

      // // Incidence
      // inc_tmc = haz_tmc * surv_lag;
      
      // Cumulative incidence
      // cum_inc_tmc = IntMat1 * inc_tmc;
    };

    // for model with no type (so only MMC)
    // TODO: Repitition here from function for MMC, redesign (with template?) to avoid this
    void calc_haz(SparseMatrix<Type> X_fixed, 
                  SparseMatrix<Type> X_time,
                  SparseMatrix<Type> X_age, 
                  SparseMatrix<Type> X_space,
                  SparseMatrix<Type> X_agetime, 
                  SparseMatrix<Type> X_agespace,
                  SparseMatrix<Type> X_spacetime, 
                  // Integration matrix
                  density::SparseMatrix<Type> IntMat1,
                  density::SparseMatrix<Type> IntMat2,
                  // parameters
                  vector<Type> u_fixed, 
                  vector<Type> u_age,
                  vector<Type> u_time,
                  vector<Type> u_space, 
                  array<Type> u_agetime,
                  array<Type> u_agespace,
                  array<Type> u_spacetime,
                  Type sigma_age,
                  Type sigma_time,
                  Type sigma_space,
                  Type sigma_agetime,
                  Type sigma_agespace,
                  Type sigma_spacetime) {

      // Vector the interaction terms
      vector<Type> u_agespace_v(u_agespace);
      vector<Type> u_agetime_v(u_agetime);
      vector<Type> u_spacetime_v(u_spacetime);

      /// Estimate hazard rate ///
      /// TODO: Break this down into functions as well!
      // Medical hazard rate
      haz = X_fixed * u_fixed +
        X_age * u_age * sigma_age +
        X_space * u_space * sigma_space +
				X_time * u_time * sigma_time +
				X_agetime * u_agetime_v * sigma_agetime +
				X_agespace * u_agespace_v * sigma_agespace +
				X_spacetime * u_spacetime_v * sigma_spacetime;

      // Rates on [0,1] scale
      haz = invlogit_vec(haz);

      // Survival probabilities (TODO: These only need to be calculated once!)
	    // vector<Type> logprob  = log(Type(1.0) - haz);
	    // vector<Type> surv     = exp(IntMat1 * logprob);
	    // vector<Type> surv_lag = exp(IntMat2 * logprob);
	    // leftcens              = Type(1.0) - surv;

      // Incidence
      // inc = haz * surv_lag;

      // Cumulative incidence
      // cum_inc = IntMat1 * inc;
    };

    // final calculation of total report vals (e.g. haz = haz_mmc + haz_tmc)
    void calc_haz() {

      // Adjustment such that \lambda_mmc + \lambda_tmc \in [0,1]  
      // Medical rate to only take from the remaining proportion 
      // not taken through traditional circumcision (1 - \lambda_tmc)
      haz_mmc = haz_mmc * (1 - haz_tmc);
      
      // Total hazard rate
      haz = haz_mmc + haz_tmc;
    };



    // final calculation of total report vals (e.g. haz = haz_mmc + haz_tmc)
    void calc_report_vals(// Integration matrices
                          density::SparseMatrix<Type> IntMat1,
                          density::SparseMatrix<Type> IntMat2) {

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
      inc = haz * surv_lag;

      // Cumulative incidence
      cum_inc = cum_inc_tmc + cum_inc_mmc;
    };

  // function to calculate survival probabilities
  void calc_surv(// Integration matrices
                  density::SparseMatrix<Type> IntMat1,
                  density::SparseMatrix<Type> IntMat2) {

      logprob  = log(Type(1.0) - haz);
	    surv     = exp(IntMat1 * logprob);
	    surv_lag = exp(IntMat2 * logprob);
	    leftcens = Type(1.0) - surv;
  }

    // Function to calculate incidence & cumulative incidence
    // TODO: This code is reused, find a way to change e.g inc_mmc while referring to inc here
    // Actually works well here with if statement, but may not scale well with more models
    void calc_inc(density::SparseMatrix<Type> IntMat1, int is_type) {

      // Incidence
      if (is_type == 1) {
        inc_mmc = haz_mmc * surv_lag;
        inc_tmc = haz_tmc * surv_lag;
      }
      inc = haz * surv_lag; // TODO: Ask Matt why this isn't inc_mmc + inc_tmc for model w/ type??

      // Cumulative incidence
      if (is_type == 1) {
        cum_inc_mmc = IntMat1 * inc_mmc;
        cum_inc_tmc = IntMat1 * inc_tmc;
        cum_inc     = cum_inc_tmc + cum_inc_mmc;
      } else {
        cum_inc = IntMat1 * inc;
      }
    }

    // Function to calculate likelihood
    // TODO: Can make an enum (or something?) pointer to iterate over for this
    // (will be different for each model)
  void likelihood(density::SparseMatrix<Type> Mat,
                  vector<Type> report_val) {
      // Calculate likelihood for chosen circ pop with corresponding report vals
    nll -= (Mat * log(report_val)).sum();
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
    vector<Type> get_logprob() {
      return logprob;
    };
    vector<Type> get_surv() {
      return surv;
    };
    vector<Type> get_surv_lag() {
      return surv_lag;
    };
    vector<Type> get_leftcens() {
      return leftcens;
    };
};

#endif

