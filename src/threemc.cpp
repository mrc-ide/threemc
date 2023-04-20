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
  // TODO: Change all of these nested ifs into a switch statement, very ugly!!
  if (threemc_data.is_type == 1) {
    if (threemc_data.rw_order == 1) {
      Threemc_rw<Type> threemc;
      // std::cout << "Message about model selection" << std::endl;
      threemc.calc_nll(threemc_data, this);
      nll = threemc.get_nll();
    } else {
      Threemc<Type> threemc;
      // std::cout << "Message about model selection" << std::endl;
      threemc.calc_nll(threemc_data, this);
      nll = threemc.get_nll();
     }
    // class for model with no type
  } else if (threemc_data.is_type == 0) {
    if (threemc_data.rw_order == 1) {
      Threemc_nt_rw<Type> threemc;
    //   // std::cout << "Message about model selection" << std::endl;
      threemc.calc_nll(threemc_data, this);
      nll = threemc.get_nll();
    } else {
      Threemc_nt<Type> threemc;
      // std::cout << "Message about model selection" << std::endl;
      threemc.calc_nll(threemc_data, this);
      nll = threemc.get_nll();
    }
  }
  return nll;
}
