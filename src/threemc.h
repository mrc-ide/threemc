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

// Data Struct //
template<class Type>
struct Threemc_data {

  // Survival analysis matrices
  density::SparseMatrix<Type> A_mmc; // Matrix selecting instantaneous hazard for medically circumcised pop
  density::SparseMatrix<Type> A_tmc; // Matrix selecting instantaneous hazard for traditionally circumcised pop
  density::SparseMatrix<Type> A_mc; // Matrix selecting instantaneous hazard for unknown circumcised pop
  density::SparseMatrix<Type> B; // Matrix selecting relevant cumulative hazard entry for observed and right censored pop
  density::SparseMatrix<Type> C; // Matrix selecting relevant cumulative hazard entry for interval censored pop

  // Integeration matrices
  density::SparseMatrix<Type> IntMat1;
  density::SparseMatrix<Type> IntMat2;

  // Design Matrices
  density::SparseMatrix<Type> X_fixed_mmc; // Design matrix for the fixed effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_time_mmc; // Design matrix for the temporal random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_age_mmc; // Design matrix for the stratification random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_space_mmc; // Design matrix for the stratification random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_agetime_mmc; // Design matrix for the interaction random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_agespace_mmc; // Design matrix for the interaction random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_spacetime_mmc; // Design matrix for the interaction random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_fixed_tmc; // Design matrix for the fixed effects in the traditional circumcision hazard rate
  density::SparseMatrix<Type> X_age_tmc; // Design matrix for the stratification random effects in the traditional circumcision hazard rate
  density::SparseMatrix<Type> X_space_tmc; // Design matrix for the stratification random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_agespace_tmc; // Design matrix for the interaction random effects in the medical circumcision hazard rate

  // for model with no type
  density::SparseMatrix<Type> X_fixed;    // Design matrix for the fixed effects
  density::SparseMatrix<Type> X_time;     // Design matrix for the temporal random effects
  density::SparseMatrix<Type> X_age;      // Design matrix for the stratification random effects
  density::SparseMatrix<Type> X_space;    // Design matrix for the stratification random effects
  density::SparseMatrix<Type> X_agetime;  // Design matrix for the age-time interaction random effects
  density::SparseMatrix<Type> X_agespace; // Design matrix for the age-space interaction random effects
  density::SparseMatrix<Type> X_spacetime; // Design matrix for the age-space interaction random effects
 
  // Precision Matrix
  density::SparseMatrix<Type> Q_space; // Aggregation matrix for number of circumcisions performed

  // Constructor, which conditionally uses data macros to assign values to declared matrices
  // Threemc_data(int is_type) {
  Threemc_data(SEXP x) {

    // common to all
    A_mc    = tmbutils::asSparseMatrix<Type>(getListElement(x, "A_mc"));
    B       = tmbutils::asSparseMatrix<Type>(getListElement(x, "B"));
    C       = tmbutils::asSparseMatrix<Type>(getListElement(x, "C"));
    IntMat1 = tmbutils::asSparseMatrix<Type>(getListElement(x, "IntMat1"));
    IntMat2 = tmbutils::asSparseMatrix<Type>(getListElement(x, "IntMat2"));

    Q_space = tmbutils::asSparseMatrix<Type>(getListElement(x, "Q_space"));

    A_mmc = tmbutils::asSparseMatrix<Type>(getListElement(x, "A_mmc"));
    A_tmc = tmbutils::asSparseMatrix<Type>(getListElement(x, "A_tmc"));

    X_fixed_mmc    = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_fixed_mmc"));
    X_time_mmc     = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_time_mmc")); 
    X_age_mmc      = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_age_mmc")); 
    X_space_mmc    = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_space_mmc")); 
    X_agetime_mmc  = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_agetime_mmc")); 
    X_agespace_mmc = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_agespace_mmc")); 
    X_spacetime_mmc= tmbutils::asSparseMatrix<Type>(getListElement(x, "X_spacetime_mmc"));
    X_fixed_tmc    = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_fixed_tmc")); 
    X_age_tmc      = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_age_tmc")); 
    X_space_tmc    = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_space_tmc")); 
    X_agespace_tmc = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_agespace_tmc"));
  }
};


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
    // Need to just overload this function
    // TODO: Can definitely design this function better to avoid repitition
    // For MMC: 
    void calc_haz(density::SparseMatrix<Type> X_fixed, 
                  density::SparseMatrix<Type> X_time,
                  density::SparseMatrix<Type> X_age, 
                  density::SparseMatrix<Type> X_space,
                  density::SparseMatrix<Type> X_agetime, 
                  density::SparseMatrix<Type> X_agespace,
                  density::SparseMatrix<Type> X_spacetime, 
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
    };

    // For TMC: 
    void calc_haz(density::SparseMatrix<Type> X_fixed,
                  density::SparseMatrix<Type> X_age, 
                  density::SparseMatrix<Type> X_space,
                  density::SparseMatrix<Type> X_agespace,
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
      // if (is_type == 1) {
        inc_mmc = haz_mmc * surv_lag;
        inc_tmc = haz_tmc * surv_lag;
      // }
      inc = haz * surv_lag; // TODO: Ask Matt why this isn't inc_mmc + inc_tmc for model w/ type??

      // Cumulative incidence
      // if (is_type == 1) {
        cum_inc_mmc = IntMat1 * inc_mmc;
        cum_inc_tmc = IntMat1 * inc_tmc;
        cum_inc     = cum_inc_tmc + cum_inc_mmc;
      // } else {
      //   cum_inc = IntMat1 * inc;
      // }
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

