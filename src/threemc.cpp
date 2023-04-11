/// @file threemc.cpp
#define TMB_LIB_INIT R_init_threemc
#include <TMB.hpp>
#include "threemc.h"
#include "implementation.cpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace density;

  // Define struct containing data matrices
  DATA_STRUCT(threemc_data, Threemc_data);
  
  // define object of class Threemc, which will store our negative log likelihood
  Threemc<Type> threemc;

  // calculate nll and report values of interest (haz, inc, cum_inc & surv prob)
  threemc.calc_nll(threemc_data, this);
  
  //// return nll ////
  return threemc.get_nll();
}
