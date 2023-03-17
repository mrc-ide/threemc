#define TMB_LIB_INIT R_init_threemc
#include "threemc.h"
// #include "implementation.cpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace density;

  // define object of class Threemc, which will store our negative log likelihood
  Threemc<Type> threemc;

  //////////////////
  /// Parameters ///
  //////////////////

  // Fixed Effects
  PARAMETER_VECTOR(u_fixed_mmc);
  PARAMETER_VECTOR(u_fixed_tmc);
  
  // /// Calculate nll ///
  
  // Apply prior on fixed effects
  threemc.fixed_effects_prior(u_fixed_mmc, u_fixed_tmc);

  /// return nll ///
  return threemc.get_nll();
}
