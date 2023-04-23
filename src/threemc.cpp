/// @file threemc.cpp
#define TMB_LIB_INIT R_init_threemc
// #include <iostream>
#include <TMB.hpp>
#include "threemc.h"
#include "implementation.cpp"

// convert conditions to a bit, concatenate these bits into an int to switch on
#define TYPE (1 << 0)
#define RW (1 << 1)
#define PAED (1 << 2)

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace density;

  // Define struct containing data matrices
  DATA_STRUCT(threemc_data, Threemc_data);
  
  // initialise nll
  Type nll = Type(0);

  // Switch statement where particular form of model is decided
  switch((threemc_data.is_type? TYPE : 0) | (threemc_data.rw_order? RW : 0) |
         (threemc_data.paed_age_cutoff? PAED : 0)) {
    case 0: // no type
      nll = nll_switch<Type, Threemc_nt<Type>>(nll, threemc_data, this);
      break;
    case TYPE: // "default" model with type information
      nll = nll_switch<Type, Threemc<Type>>(nll, threemc_data, this);
      break;      
    case RW: // no type, RW temporal prior
      nll = nll_switch<Type, Threemc_nt_rw<Type>>(nll, threemc_data, this);
      break;      
    case TYPE + RW: // type info, RW temporal prior
      nll = nll_switch<Type, Threemc_rw<Type>>(nll, threemc_data, this);
      break;      
    case PAED: // No model for paediatric age cutoff with no type, fall through
    case TYPE + PAED: // type info, paediatric age cutoff for MMC
      nll = nll_switch<Type, Threemc_paed<Type>>(nll, threemc_data, this);
      break;
    // case RW + PAED: 
    // case TYPE + RW + PAED: // RW temporal prior, paediatric age cutoff for MMC
    //   nll = nll_switch<Type, Threemc_rw_paed<Type>>(nll, threemc_data, this);
    //   break;
  }

  return nll;
}
