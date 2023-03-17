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
    };

    // TODO: Write constructor with more arguments?

    // Default Destructor (needed??)
    // ~Threemc() {
    // };

    // Prior on fixed effects (won't need second version until implementing no type version)
    void fixed_effects_prior(vector<Type> u_fixed_mmc, vector<Type> u_fixed_tmc) {
      // fixed effects for the medical circumcision rate
      nll -= dnorm(u_fixed_mmc, Type(0), Type(5), true).sum();
      // fixed effects for the traditional circumcision rate
      nll -= dnorm(u_fixed_tmc, Type(0), Type(5), true).sum();
    };

    // getter for nll;
    Type get_nll() {
      return nll;
    };
};

#endif

