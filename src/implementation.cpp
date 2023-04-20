/// @file implementation.cpp
// #include <TMB.hpp>
// #include "threemc.h"

/* Utility Functions (perhaps move to the end? Or their own header?) */

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

//// Data Struct ////

// Constructor, which conditionally uses data macros to assign values to declared matrices
// Threemc_data(int is_type) {
template<class Type>
Threemc_data<Type>::Threemc_data(SEXP x) {

  // indicators
  // Does model include type or not? 
  is_type  = CppAD::Integer(asVector<Type>(getListElement(x, "is_type"))[0]);
  // Use AR1 (0) or RW (1) temporal prior
  rw_order = CppAD::Integer(asVector<Type>(getListElement(x, "rw_order"))[0]);
  
  // common to all
  A_mc    = tmbutils::asSparseMatrix<Type>(getListElement(x, "A_mc"));
  B       = tmbutils::asSparseMatrix<Type>(getListElement(x, "B"));
  C       = tmbutils::asSparseMatrix<Type>(getListElement(x, "C"));
  IntMat1 = tmbutils::asSparseMatrix<Type>(getListElement(x, "IntMat1"));
  IntMat2 = tmbutils::asSparseMatrix<Type>(getListElement(x, "IntMat2"));

  Q_space = tmbutils::asSparseMatrix<Type>(getListElement(x, "Q_space"));

  // for model with type only
  if (is_type == 1) {
    A_mmc = tmbutils::asSparseMatrix<Type>(getListElement(x, "A_mmc"));
    A_tmc = tmbutils::asSparseMatrix<Type>(getListElement(x, "A_tmc"));

    X_fixed_mmc     = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_fixed_mmc"));
    X_time_mmc      = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_time_mmc")); 
    X_age_mmc       = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_age_mmc")); 
    X_space_mmc     = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_space_mmc")); 
    X_agetime_mmc   = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_agetime_mmc")); 
    X_agespace_mmc  = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_agespace_mmc")); 
    X_spacetime_mmc = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_spacetime_mmc"));
    X_fixed_tmc     = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_fixed_tmc")); 
    X_age_tmc       = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_age_tmc")); 
    X_space_tmc     = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_space_tmc")); 
    X_agespace_tmc  = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_agespace_tmc")); 

  // for model with no type
  } else {
    X_fixed     = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_fixed"));
    X_time      = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_time")); 
    X_age       = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_age")); 
    X_space     = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_space")); 
    X_agetime   = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_agetime")); 
    X_agespace  = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_agespace")); 
    X_spacetime = tmbutils::asSparseMatrix<Type>(getListElement(x, "X_spacetime"));
  }

  if (rw_order == 1) {
    Q_time = tmbutils::asSparseMatrix<Type>(getListElement(x, "Q_time"));
  }
}


//// Threemc class ////

// Constructor
template<class Type> 
Threemc<Type>::Threemc() {
  Type nll = Type(0);
}

// Destructor
template<class Type> 
Threemc<Type>::~Threemc() {
}

// Fixed effects for circumcision rate
template<class Type> 
void Threemc<Type>::fix_eff_p(vector<Type> u_fixed) {
  nll -= dnorm(u_fixed, Type(0), Type(5), true).sum();
}

// Prior on temporal random effects
// This function works for both type and no type specifications (just with different inputs)
template<class Type> 
void Threemc<Type>::rand_eff_time_p(vector<Type> u_time,
                                    Type logsigma_time,
                                    Type sigma_time,
                                    Type logitrho_time1,
                                    Type rho_time1) {

  // density::AR1 Process
  nll += density::AR1(rho_time1)(u_time);

  // Sum to zero constraint
  nll -= dnorm(u_time.sum(), Type(0), Type(0.001) * u_time.size(), true);

  // Prior on the standard deviation for the temporal random effects
  nll -= dexp(sigma_time, Type(1), true) + logsigma_time;

  // Prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_time1, Type(3), Type(3), true);
}

// Prior on the age random effects
// TODO: Will change in model where there is a paediatric-adult age split
template<class Type>
void Threemc<Type>::rand_eff_age_p(vector<Type> u_age,
                                   Type logsigma_age,
                                   Type sigma_age,
                                   Type logitrho_age1,
                                   Type rho_age1) {

  // density::AR1 processes
  nll += density::AR1(rho_age1)(u_age);

  // sum to zero constraint
  nll -= dnorm(u_age.sum(), Type(0), Type(0.001) * u_age.size(), true);

  // prior on the standard deviation for the age random effects
  nll -= dexp(sigma_age, Type(1), true) + logsigma_age;

  // prior on the logit autocorrelation parameters
  nll -= dnorm(logitrho_age1, Type(3), Type(2), true);
}

// Prior on the spatial random effects
template<class Type>
void Threemc<Type>::rand_eff_space_p(density::SparseMatrix<Type> Q_space,
                                     vector<Type> u_space,
                                     Type logsigma_space,
                                     Type sigma_space) {

  // Gaussian markov random field with prespecified precision matrix
  nll += GMRF(Q_space)(u_space);

  // Sum to zero constraints
  nll -= dnorm(u_space.sum(), Type(0), Type(0.001) * u_space.size(), true);

  // Prior on the standard deviation for the spatial random effects
  nll -= dexp(sigma_space, Type(1), true) + logsigma_space;
}

// Prior on the interaction random effects for either MMC or MC (for model w/ no type)
template<class Type>
void Threemc<Type>::rand_eff_interact_p(density::SparseMatrix<Type> Q_space,
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

  // Interactions: space-time (GMRF x density::AR1), age-time (density::AR1 x density::AR1) and age-space (density::AR1 x GMRF)
  nll += SEPARABLE(density::AR1(rho_time2), density::AR1(rho_age2))(u_agetime);
  nll += SEPARABLE(GMRF(Q_space), density::AR1(rho_age3))(u_agespace);
  nll += SEPARABLE(GMRF(Q_space), density::AR1(rho_time3))(u_spacetime);
  
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
}

// Overload prior on the interaction random effects for just TMC
template<class Type>
void Threemc<Type>::rand_eff_interact_p(density::SparseMatrix<Type> Q_space,
                                        array<Type> u_agespace_tmc,
                                        Type logsigma_agespace_tmc,
                                        Type sigma_agespace_tmc,
                                        Type logitrho_tmc_age2,
                                        Type rho_tmc_age2) {

  // Interactions: space-time (GMRF x density::AR1), age-time (density::AR1 x density::AR1) and age-space (density::AR1 x GMRF)
  nll += SEPARABLE(GMRF(Q_space), density::AR1(rho_tmc_age2))(u_agespace_tmc);
  
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
}

// Function to calculate report values 
// Need to just overload this function
// TODO: Can definitely design this function better to avoid repitition
// For MMC: 
template<class Type> 
void Threemc<Type>::calc_haz(density::SparseMatrix<Type> X_fixed, 
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
}

// For TMC: 
template<class Type> 
void Threemc<Type>::calc_haz(density::SparseMatrix<Type> X_fixed, 
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
}

// final calculation of total report vals (e.g. haz = haz_mmc + haz_tmc)
template<class Type>
void Threemc<Type>::calc_haz() {

  // Adjustment such that \lambda_mmc + \lambda_tmc \in [0,1]  
  // Medical rate to only take from the remaining proportion 
  // not taken through traditional circumcision (1 - \lambda_tmc)
  haz_mmc = haz_mmc * (1 - haz_tmc);
  
  // Total hazard rate
  haz = haz_mmc + haz_tmc;
}

// function to calculate survival probabilities
template<class Type>
void Threemc<Type>::calc_surv(// Integration matrices
                              density::SparseMatrix<Type> IntMat1,
                              density::SparseMatrix<Type> IntMat2) {

    vector<Type> logprob  = log(Type(1.0) - haz);
    surv     = exp(IntMat1 * logprob);
    surv_lag = exp(IntMat2 * logprob);
    leftcens = Type(1.0) - surv;
}

// Function to calculate incidence & cumulative incidence
// TODO: This code is reused, find a way to change e.g inc_mmc while referring to inc here
// Actually works well here with if statement, but may not scale well with more models
template<class Type>
void Threemc<Type>::calc_inc(density::SparseMatrix<Type> IntMat1, int is_type) {
// void Threemc<Type>::calc_inc(density::SparseMatrix<Type> IntMat1) {

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

// template<class Type>
// void Threemc<Type>::calc_inc(density::SparseMatrix<Type> IntMat1, int is_type) {
//   inc = haz * surv_lag;
//   if (is_type == 1) {
//     cum_inc = IntMat1 * inc;
//   }
// }

// Function to calculate likelihood
// TODO: Can make an enum (or something?) pointer to iterate over for this
// (will be different for each model)
template<class Type> 
void Threemc<Type>::likelihood(density::SparseMatrix<Type> Mat,
                               vector<Type> report_val) {
  // Calculate likelihood for chosen circ pop with corresponding report vals
  nll -= (Mat * log(report_val)).sum();
}

//// Threem_nt class, Threemc with no type information ////

// Constructor
template<class Type> 
Threemc_nt<Type>::Threemc_nt() {
  Type nll = Type(0);
}

// Destructor
template<class Type> 
Threemc_nt<Type>::~Threemc_nt() {
}

// Calculate haz for model with no type (so only MMC)
// TODO: Repitition here from function for MMC, redesign (with template?) to avoid this
template<class Type>
void Threemc_nt<Type>::calc_haz(density::SparseMatrix<Type> X_fixed, 
                                density::SparseMatrix<Type> X_time,
                                density::SparseMatrix<Type> X_age, 
                                density::SparseMatrix<Type> X_space,
                                density::SparseMatrix<Type> X_agetime, 
                                density::SparseMatrix<Type> X_agespace,
                                density::SparseMatrix<Type> X_spacetime, 
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
}


//// Threemc_rw class, Threemc with RW temporal prior ////

// Constructor
template<class Type> 
Threemc_rw<Type>::Threemc_rw() {
  Type nll = Type(0);
}

// Destructor
template<class Type> 
Threemc_rw<Type>::~Threemc_rw() {
}

// Prior on temporal random effects
// This function works for both type and no type specifications (just with different inputs)
template<class Type> 
void Threemc_rw<Type>::rand_eff_time_p(density::SparseMatrix<Type> Q_time,
                                       vector<Type> u_time,
                                       Type logsigma_time,
                                       Type sigma_time) {
                                       // Type logitrho_time1,
                                       // Type rho_time1) {

  // density::AR1 Process
  // nll += density::AR1(rho_time1)(u_time);
  nll += GMRF(Q_time)(u_time);

  // Sum to zero constraint
  nll -= dnorm(u_time.sum(), Type(0), Type(0.001) * u_time.size(), true);

  // Prior on the standard deviation for the temporal random effects
  nll -= dexp(sigma_time, Type(1), true) + logsigma_time;

  // Prior on the logit autocorrelation parameters
  // nll -= dnorm(logitrho_time1, Type(3), Type(3), true);
}

// Prior on the interaction random effects for either MMC or MC (for model w/ no type)
template<class Type>
void Threemc_rw<Type>::rand_eff_interact_p(density::SparseMatrix<Type> Q_space,
                                           density::SparseMatrix<Type> Q_time,
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
                                           Type rho_age3) {
                                           // Type logitrho_time2,
                                           // Type rho_time2,
                                           // Type logitrho_time3,
                                           // Type rho_time3) {

  
  // Interactions: space-time (GMRF x density::AR1), age-time (density::AR1 x density::AR1) and age-space (density::AR1 x GMRF)
  nll += SEPARABLE(GMRF(Q_time), density::AR1(rho_age2))(u_agetime);
  nll += SEPARABLE(GMRF(Q_space), density::AR1(rho_age3))(u_agespace);
  nll += SEPARABLE(GMRF(Q_space), GMRF(Q_time))(u_spacetime);
  
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
  nll -= dnorm(logitrho_age2,  Type(3), Type(2), true);
  nll -= dnorm(logitrho_age3,  Type(3), Type(2), true);
}

//// Threemc_nw_rw class, Threemc with no type split & RW temporal prior ////

// Constructor
template<class Type> 
Threemc_nt_rw<Type>::Threemc_nt_rw() {
  Type nll = Type(0);
}

// Destructor
template<class Type> 
Threemc_nt_rw<Type>::~Threemc_nt_rw() {
}


//// Functions to pull through parameters & calculate NLL ////
  
// redefine TMB_OBJECTIVE_PTR so that members can access parameters
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// Function to calculate NLL and report values of interest
template<class Type> 
void Threemc<Type>::calc_nll(struct Threemc_data<Type> threemc_data,
                             objective_function<Type>* obj) {
  // Parameters

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
  PARAMETER(logsigma_age_mmc);       
  PARAMETER(logsigma_time_mmc);      
  PARAMETER(logsigma_space_mmc);     
  PARAMETER(logsigma_agetime_mmc);   
  PARAMETER(logsigma_agespace_mmc);  
  PARAMETER(logsigma_spacetime_mmc); 
  PARAMETER(logsigma_age_tmc);       
  PARAMETER(logsigma_space_tmc);     
  PARAMETER(logsigma_agespace_tmc);  

  Type sigma_age_mmc       = exp(logsigma_age_mmc);
  Type sigma_time_mmc      = exp(logsigma_time_mmc);
  Type sigma_space_mmc     = exp(logsigma_space_mmc);
  Type sigma_agetime_mmc   = exp(logsigma_agetime_mmc);
  Type sigma_agespace_mmc  = exp(logsigma_agespace_mmc);
  Type sigma_spacetime_mmc = exp(logsigma_spacetime_mmc);
  Type sigma_age_tmc       = exp(logsigma_age_tmc);
  Type sigma_space_tmc     = exp(logsigma_space_tmc);
  Type sigma_agespace_tmc  = exp(logsigma_agespace_tmc);

  // Autocorrelation parameters 
  PARAMETER(logitrho_mmc_time1);
  PARAMETER(logitrho_mmc_time2);
  PARAMETER(logitrho_mmc_time3);
  PARAMETER(logitrho_mmc_age1);
  PARAMETER(logitrho_mmc_age2);
  PARAMETER(logitrho_mmc_age3);
  PARAMETER(logitrho_tmc_age1);
  PARAMETER(logitrho_tmc_age2);

  Type rho_mmc_time1  = geninvlogit(logitrho_mmc_time1, Type(-1.0), Type(1.0));
  Type rho_mmc_time2  = geninvlogit(logitrho_mmc_time2, Type(-1.0), Type(1.0));
  Type rho_mmc_time3  = geninvlogit(logitrho_mmc_time3, Type(-1.0), Type(1.0));
  Type rho_mmc_age1   = geninvlogit(logitrho_mmc_age1,  Type(-1.0), Type(1.0));
  Type rho_mmc_age2   = geninvlogit(logitrho_mmc_age2,  Type(-1.0), Type(1.0));
  Type rho_mmc_age3   = geninvlogit(logitrho_mmc_age3,  Type(-1.0), Type(1.0));
  Type rho_tmc_age1   = geninvlogit(logitrho_tmc_age1,  Type(-1.0), Type(1.0));
  Type rho_tmc_age2   = geninvlogit(logitrho_tmc_age2,  Type(-1.0), Type(1.0));

  //// Priors ////

  // Apply prior on fixed effects
  fix_eff_p(u_fixed_mmc);
  fix_eff_p(u_fixed_tmc);

  // prior on temporal random effects
  rand_eff_time_p(u_time_mmc,
                  logsigma_time_mmc,
                  sigma_time_mmc,
                  logitrho_mmc_time1,
                  rho_mmc_time1);

  // prior on the age random effects
  rand_eff_age_p(u_age_mmc,
                 logsigma_age_mmc,
                 sigma_age_mmc,
                 logitrho_mmc_age1,
                 rho_mmc_age1);
  rand_eff_age_p(u_age_tmc,
                 logsigma_age_tmc,
                 sigma_age_tmc,
                 logitrho_tmc_age1,
                 rho_tmc_age1);

  // prior on spatial random effects
  rand_eff_space_p(threemc_data.Q_space,
                   u_space_mmc,
                   logsigma_space_mmc,
                   sigma_space_mmc);
  rand_eff_space_p(threemc_data.Q_space,
                   u_space_tmc,
                   logsigma_space_tmc,
                   sigma_space_tmc);

  // prior on interaction random effects
  rand_eff_interact_p(threemc_data.Q_space,
                      u_agespace_mmc,
                      u_agetime_mmc,
                      u_spacetime_mmc,
                      logsigma_agespace_mmc,
                      sigma_agespace_mmc,
                      logsigma_agetime_mmc,
                      sigma_agetime_mmc,
                      logsigma_spacetime_mmc,
                      sigma_spacetime_mmc,
                      logitrho_mmc_age2,
                      rho_mmc_age2,
                      logitrho_mmc_age3,
                      rho_mmc_age3,
                      logitrho_mmc_time2,
                      rho_mmc_time2,
                      logitrho_mmc_time3,
                      rho_mmc_time3);
  rand_eff_interact_p(threemc_data.Q_space,
                      u_agespace_tmc,
                      logsigma_agespace_tmc,
                      sigma_agespace_tmc,
                      logitrho_tmc_age2,
                      rho_tmc_age2);

  //// Calculate report values (hazard, (cumulative) incidence) ////

  // Calculate hazards
  calc_haz(threemc_data.X_fixed_mmc, 
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

  calc_haz(threemc_data.X_fixed_tmc, 
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

  calc_haz();

  // calculate survival probabilities
  calc_surv(threemc_data.IntMat1, threemc_data.IntMat2);

  // calculate incidences and cumulative incidences
  calc_inc(threemc_data.IntMat1, threemc_data.is_type); 

  //// Calculate likelihood ////
  likelihood(threemc_data.A_mmc, inc_mmc); // medical circumcisions
  likelihood(threemc_data.A_tmc, inc_tmc); // traditional circumcisions
  likelihood(threemc_data.A_mc, inc);      // circs of unknown type
  likelihood(threemc_data.B, surv);        // right censored (i.e. uncircumcised)
  likelihood(threemc_data.C, leftcens);    // left censored (i.e. unknown circ age)

  //// report hazard rates, incidence and cumulative incidence ////
  REPORT(haz_mmc);     // Medical hazard rate
  REPORT(haz_tmc);     // Traditional hazard rate
  REPORT(haz);         // Total hazard rate
  REPORT(inc_mmc);     // Medical circumcision incidence rate
  REPORT(inc_tmc);     // Traditional circumcision incidence rate
  REPORT(inc);         // Total circumcision incidence rate
  REPORT(cum_inc_mmc); // Medical circumcision cumulative incidence rate
  REPORT(cum_inc_tmc); // Traditional circumcision cumulative incidence rate
  REPORT(cum_inc);     // Total circumcision cumulative incidence rate
  REPORT(surv);        // Survival probabilities
}

// For model with RW temporal prior
template<class Type> 
void Threemc_rw<Type>::calc_nll(struct Threemc_data<Type> threemc_data,
                             objective_function<Type>* obj) {

  // Parameters

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
  PARAMETER(logsigma_age_mmc);       
  PARAMETER(logsigma_time_mmc);      
  PARAMETER(logsigma_space_mmc);     
  PARAMETER(logsigma_agetime_mmc);   
  PARAMETER(logsigma_agespace_mmc);  
  PARAMETER(logsigma_spacetime_mmc); 
  PARAMETER(logsigma_age_tmc);       
  PARAMETER(logsigma_space_tmc);     
  PARAMETER(logsigma_agespace_tmc);  

  Type sigma_age_mmc       = exp(logsigma_age_mmc);
  Type sigma_time_mmc      = exp(logsigma_time_mmc);
  Type sigma_space_mmc     = exp(logsigma_space_mmc);
  Type sigma_agetime_mmc   = exp(logsigma_agetime_mmc);
  Type sigma_agespace_mmc  = exp(logsigma_agespace_mmc);
  Type sigma_spacetime_mmc = exp(logsigma_spacetime_mmc);
  Type sigma_age_tmc       = exp(logsigma_age_tmc);
  Type sigma_space_tmc     = exp(logsigma_space_tmc);
  Type sigma_agespace_tmc  = exp(logsigma_agespace_tmc);

  // Autocorrelation parameters 
  // PARAMETER(logitrho_mmc_time1);
  // PARAMETER(logitrho_mmc_time2);
  // PARAMETER(logitrho_mmc_time3);
  PARAMETER(logitrho_mmc_age1);
  PARAMETER(logitrho_mmc_age2);
  PARAMETER(logitrho_mmc_age3);
  PARAMETER(logitrho_tmc_age1);
  PARAMETER(logitrho_tmc_age2);

  // Type rho_mmc_time1  = geninvlogit(logitrho_mmc_time1, Type(-1.0), Type(1.0));
  // Type rho_mmc_time2  = geninvlogit(logitrho_mmc_time2, Type(-1.0), Type(1.0));
  // Type rho_mmc_time3  = geninvlogit(logitrho_mmc_time3, Type(-1.0), Type(1.0));
  Type rho_mmc_age1   = geninvlogit(logitrho_mmc_age1,  Type(-1.0), Type(1.0));
  Type rho_mmc_age2   = geninvlogit(logitrho_mmc_age2,  Type(-1.0), Type(1.0));
  Type rho_mmc_age3   = geninvlogit(logitrho_mmc_age3,  Type(-1.0), Type(1.0));
  Type rho_tmc_age1   = geninvlogit(logitrho_tmc_age1,  Type(-1.0), Type(1.0));
  Type rho_tmc_age2   = geninvlogit(logitrho_tmc_age2,  Type(-1.0), Type(1.0));

  //// Priors ////

  // Apply prior on fixed effects
  fix_eff_p(u_fixed_mmc);
  fix_eff_p(u_fixed_tmc);

  // prior on temporal random effects
  rand_eff_time_p(threemc_data.Q_time,
                  u_time_mmc,
                  logsigma_time_mmc,
                  sigma_time_mmc);

  // prior on the age random effects
  rand_eff_age_p(u_age_mmc,
                 logsigma_age_mmc,
                 sigma_age_mmc,
                 logitrho_mmc_age1,
                 rho_mmc_age1);
  rand_eff_age_p(u_age_tmc,
                 logsigma_age_tmc,
                 sigma_age_tmc,
                 logitrho_tmc_age1,
                 rho_tmc_age1);

  // prior on spatial random effects
  rand_eff_space_p(threemc_data.Q_space,
                   u_space_mmc,
                   logsigma_space_mmc,
                   sigma_space_mmc);
  rand_eff_space_p(threemc_data.Q_space,
                   u_space_tmc,
                   logsigma_space_tmc,
                   sigma_space_tmc);

  // prior on interaction random effects
  rand_eff_interact_p(threemc_data.Q_space,
                      threemc_data.Q_time,
                      u_agespace_mmc,
                      u_agetime_mmc,
                      u_spacetime_mmc,
                      logsigma_agespace_mmc,
                      sigma_agespace_mmc,
                      logsigma_agetime_mmc,
                      sigma_agetime_mmc,
                      logsigma_spacetime_mmc,
                      sigma_spacetime_mmc,
                      logitrho_mmc_age2,
                      rho_mmc_age2,
                      logitrho_mmc_age3,
                      rho_mmc_age3);

  rand_eff_interact_p(threemc_data.Q_space,
                      u_agespace_tmc,
                      logsigma_agespace_tmc,
                      sigma_agespace_tmc,
                      logitrho_tmc_age2,
                      rho_tmc_age2);

  //// Calculate report values (hazard, (cumulative) incidence) ////

  // Calculate hazards
  calc_haz(threemc_data.X_fixed_mmc, 
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

  calc_haz(threemc_data.X_fixed_tmc, 
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

  calc_haz();

  // calculate survival probabilities
  calc_surv(threemc_data.IntMat1, threemc_data.IntMat2);

  // calculate incidences and cumulative incidences
  calc_inc(threemc_data.IntMat1, threemc_data.is_type); 

  //// Calculate likelihood ////
  likelihood(threemc_data.A_mmc, inc_mmc); // medical circumcisions
  likelihood(threemc_data.A_tmc, inc_tmc); // traditional circumcisions
  likelihood(threemc_data.A_mc, inc);      // circs of unknown type
  likelihood(threemc_data.B, surv);        // right censored (i.e. uncircumcised)
  likelihood(threemc_data.C, leftcens);    // left censored (i.e. unknown circ age)

  //// report hazard rates, incidence and cumulative incidence ////
  REPORT(haz_mmc);     // Medical hazard rate
  REPORT(haz_tmc);     // Traditional hazard rate
  REPORT(haz);         // Total hazard rate
  REPORT(inc_mmc);     // Medical circumcision incidence rate
  REPORT(inc_tmc);     // Traditional circumcision incidence rate
  REPORT(inc);         // Total circumcision incidence rate
  REPORT(cum_inc_mmc); // Medical circumcision cumulative incidence rate
  REPORT(cum_inc_tmc); // Traditional circumcision cumulative incidence rate
  REPORT(cum_inc);     // Total circumcision cumulative incidence rate
  REPORT(surv);        // Survival probabilities
}

template<class Type>
void Threemc_nt<Type>::calc_nll(struct Threemc_data<Type> threemc_data,
                                objective_function<Type>* obj) {

  // Parameters
    
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

  Type sigma_age       = exp(logsigma_age);
  Type sigma_time      = exp(logsigma_time);
  Type sigma_space     = exp(logsigma_space);
  Type sigma_agetime   = exp(logsigma_agetime);
  Type sigma_agespace  = exp(logsigma_agespace);
  Type sigma_spacetime = exp(logsigma_spacetime);

  // Autocorrelation parameters 
  PARAMETER(logitrho_time1);
  PARAMETER(logitrho_time2);
  PARAMETER(logitrho_time3);
  PARAMETER(logitrho_age1);
  PARAMETER(logitrho_age2);
  PARAMETER(logitrho_age3);

  Type rho_time1  = geninvlogit(logitrho_time1, Type(-1.0), Type(1.0));
  Type rho_time2  = geninvlogit(logitrho_time2, Type(-1.0), Type(1.0));
  Type rho_time3  = geninvlogit(logitrho_time3, Type(-1.0), Type(1.0));
  Type rho_age1   = geninvlogit(logitrho_age1,  Type(-1.0), Type(1.0));
  Type rho_age2   = geninvlogit(logitrho_age2,  Type(-1.0), Type(1.0));
  Type rho_age3   = geninvlogit(logitrho_age3,  Type(-1.0), Type(1.0));

  //// Priors ////

  // Apply prior on fixed effects
  fix_eff_p(u_fixed);

  // prior on temporal random effects
  rand_eff_time_p(u_time,
                  logsigma_time,
                  sigma_time,
                  logitrho_time1,
                  rho_time1);

  // prior on the age random effects
  rand_eff_age_p(u_age,
                 logsigma_age,
                 sigma_age,
                 logitrho_age1,
                 rho_age1);
   
  // prior on spatial random effects
  rand_eff_space_p(threemc_data.Q_space,
                   u_space,
                   logsigma_space,
                   sigma_space);
    
  // prior on interaction random effects
  rand_eff_interact_p(threemc_data.Q_space,
                      u_agespace,
                      u_agetime,
                      u_spacetime,
                      logsigma_agespace,
                      sigma_agespace,
                      logsigma_agetime,
                      sigma_agetime,
                      logsigma_spacetime,
                      sigma_spacetime,
                      logitrho_age2,
                      rho_age2,
                      logitrho_age3,
                      rho_age3,
                      logitrho_time2,
                      rho_time2,
                      logitrho_time3,
                      rho_time3);
  
  //// Calculate report values (hazard, (cumulative) incidence) ////
    
  // Calculate hazards
  calc_haz(threemc_data.X_fixed, 
           threemc_data.X_time,
           threemc_data.X_age, 
           threemc_data.X_space,
           threemc_data.X_agetime, 
           threemc_data.X_agespace,
           threemc_data.X_spacetime, 
           threemc_data.IntMat1,
           threemc_data.IntMat2,
           u_fixed, 
           u_age,
           u_time,
           u_space, 
           u_agetime,
           u_agespace,
           u_spacetime,
           sigma_age,
           sigma_time,
           sigma_space,
           sigma_agetime,
           sigma_agespace,
           sigma_spacetime);

  // calculate survival probabilities
  calc_surv(threemc_data.IntMat1, threemc_data.IntMat2);
  
  // calculate incidences and cumulative incidences
  calc_inc(threemc_data.IntMat1, threemc_data.is_type);
   
  //// Calculate likelihood ////
  likelihood(threemc_data.A_mc, inc);    // circs of unknown type
  likelihood(threemc_data.B, surv);      // right censored (i.e. uncircumcised)
  likelihood(threemc_data.C, leftcens);  // left censored (i.e. unknown circ age)
    
  //// report hazard rates, incidence and cumulative incidence ////
  REPORT(haz);         // Total hazard rate
  REPORT(inc);         // Total circumcision incidence rate
  REPORT(cum_inc);     // Total circumcision cumulative incidence rate
  REPORT(surv);        // Survival probabilities
}

template<class Type>
void Threemc_nt_rw<Type>::calc_nll(struct Threemc_data<Type> threemc_data,
                                   objective_function<Type>* obj) {

  // Parameters
    
  // Fixed Effects
  PARAMETER_VECTOR(u_fixed);
    
  // // Age random effect
  PARAMETER_VECTOR(u_age); 
 
  // Temporal random effects 
  PARAMETER_VECTOR(u_time);

  // Spatial random effects
  PARAMETER_VECTOR(u_space);

  // Interactions
  PARAMETER_ARRAY(u_agetime);
  PARAMETER_ARRAY(u_agespace);
  PARAMETER_ARRAY(u_spacetime);

  // // Standard deviations 
  PARAMETER(logsigma_age);       
  PARAMETER(logsigma_time);      
  PARAMETER(logsigma_space);     
  PARAMETER(logsigma_agetime);   
  PARAMETER(logsigma_agespace);  
  PARAMETER(logsigma_spacetime); 

  Type sigma_age       = exp(logsigma_age);
  Type sigma_time      = exp(logsigma_time);
  Type sigma_space     = exp(logsigma_space);
  Type sigma_agetime   = exp(logsigma_agetime);
  Type sigma_agespace  = exp(logsigma_agespace);
  Type sigma_spacetime = exp(logsigma_spacetime);

  // Autocorrelation parameters 
  PARAMETER(logitrho_age1);
  PARAMETER(logitrho_age2);
  PARAMETER(logitrho_age3);

  Type rho_age1 = geninvlogit(logitrho_age1,  Type(-1.0), Type(1.0));
  Type rho_age2 = geninvlogit(logitrho_age2,  Type(-1.0), Type(1.0));
  Type rho_age3 = geninvlogit(logitrho_age3,  Type(-1.0), Type(1.0));

  //// Priors ////

  // Apply prior on fixed effects
  fix_eff_p(u_fixed);

  // prior on temporal random effects
  rand_eff_time_p(threemc_data.Q_time,
                  u_time,
                  logsigma_time,
                  sigma_time);

  // prior on the age random effects
  rand_eff_age_p(u_age,
                 logsigma_age,
                 sigma_age,
                 logitrho_age1,
                 rho_age1);
  
  // prior on spatial random effects
  rand_eff_space_p(threemc_data.Q_space,
                   u_space,
                   logsigma_space,
                   sigma_space);

  // prior on interaction random effects
  rand_eff_interact_p(threemc_data.Q_space,
                      threemc_data.Q_time,
                      u_agespace,
                      u_agetime,
                      u_spacetime,
                      logsigma_agespace,
                      sigma_agespace,
                      logsigma_agetime,
                      sigma_agetime,
                      logsigma_spacetime,
                      sigma_spacetime,
                      logitrho_age2,
                      rho_age2,
                      logitrho_age3,
                      rho_age3);

  //// Calculate report values (hazard, (cumulative) incidence) ////
   
  // Calculate hazards
  calc_haz(threemc_data.X_fixed, 
           threemc_data.X_time,
           threemc_data.X_age, 
           threemc_data.X_space,
           threemc_data.X_agetime, 
           threemc_data.X_agespace,
           threemc_data.X_spacetime, 
           threemc_data.IntMat1,
           threemc_data.IntMat2,
           u_fixed, 
           u_age,
           u_time,
           u_space, 
           u_agetime,
           u_agespace,
           u_spacetime,
           sigma_age,
           sigma_time,
           sigma_space,
           sigma_agetime,
           sigma_agespace,
           sigma_spacetime);

  // calculate survival probabilities
  calc_surv(threemc_data.IntMat1, threemc_data.IntMat2);
  
  // calculate incidences and cumulative incidences
  calc_inc(threemc_data.IntMat1, threemc_data.is_type);
   
  //// Calculate likelihood ////
  likelihood(threemc_data.A_mc, inc);    // circs of unknown type
  likelihood(threemc_data.B, surv);      // right censored (i.e. uncircumcised)
  likelihood(threemc_data.C, leftcens);  // left censored (i.e. unknown circ age)
   
  // //// report hazard rates, incidence and cumulative incidence ////
  REPORT(haz);         // Total hazard rate
  REPORT(inc);         // Total circumcision incidence rate
  REPORT(cum_inc);     // Total circumcision cumulative incidence rate
  REPORT(surv);        // Survival probabilities
}

// Function to perform boilerplate within switch in objective fun
template<class Type, class Threemc_class>
Type test(Type nll, Threemc_data<Type> threemc_data, objective_function<Type>* obj) {
  Threemc_class threemc;
  threemc.calc_nll(threemc_data, obj);
  nll = threemc.get_nll();
  return(nll);
}

// redefine TMB_OBJECTIVE_PTR so parameters can be accessed in main function
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

