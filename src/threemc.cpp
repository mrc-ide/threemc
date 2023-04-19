/// @file threemc.cpp
#define TMB_LIB_INIT R_init_threemc
// #include <iostream>
#include <TMB.hpp>
#include "threemc.h"
#include "implementation.cpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace density;

  // Define struct containing data matrices
  DATA_STRUCT(threemc_data, Threemc_data);
  
  // initialise nll
  Type nll = Type(0);

  // model with type info, no paed age cutoff and no time TMC effect
  if (threemc_data.is_type == 1) {
    Threemc<Type> threemc;
    // std::cout << "Message about model selection" << std::endl;
    threemc.calc_nll(threemc_data, this);
    nll = threemc.get_nll();
    // class for model with no type
  } else if (threemc_data.is_type == 0) {
    Threemc_nt<Type> threemc_nt;
    // std::cout << "Message about model selection" << std::endl;
    threemc_nt.calc_nll(threemc_data, this);
    nll = threemc_nt.get_nll();
  }
  return nll;
}
