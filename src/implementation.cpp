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

//// Threemc class ////

// // Constructor
template<class Type> 
Threemc<Type>::Threemc() {
  Type nll = Type(0);
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
};






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

  // // Autocorrelation parameters 
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

  // // // prior on the age random effects
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

  // // prior on spatial random effects
  rand_eff_space_p(threemc_data.Q_space,
                   u_space_mmc,
                   logsigma_space_mmc,
                   sigma_space_mmc);
  rand_eff_space_p(threemc_data.Q_space,
                   u_space_tmc,
                   logsigma_space_tmc,
                   sigma_space_tmc);

  // // prior on interaction random effects
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
  calc_inc(threemc_data.IntMat1, 1); // run where we have type

  // // //// Calculate likelihood ////
  likelihood(threemc_data.A_mmc, inc_mmc); // medical circumcisions
  likelihood(threemc_data.A_tmc, inc_tmc); // traditional circumcisions
  likelihood(threemc_data.A_mc, inc);      // circs of unknown type
  likelihood(threemc_data.B, surv);        // right censored (i.e. uncircumcised)
  likelihood(threemc_data.C, leftcens);    // left censored (i.e. unknown circ age)

  // //// report hazard rates, incidence and cumulative incidence ////
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

// redefine TMB_OBJECTIVE_PTR so parameters can be accessed in main function
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

