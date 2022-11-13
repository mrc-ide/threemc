/// @file misc.h

/// @file util.h
/* Odd Functions, Including Report Structure */

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
/* Struct to contain items to Report                               */
/*******************************************************************/
template <class Type>
struct report_values {
  
  // Report vectors 
  vector<Type> u_fixed_mmc; // change this to be indicator of model to use
  vector<Type> haz_mmc;
  vector<Type> haz_tmc;
  vector<Type> haz;
  vector<Type> inc_mmc;
  vector<Type> inc_tmc;
  vector<Type> inc;
  vector<Type> cum_inc_mmc;
  vector<Type> cum_inc_tmc;
  vector<Type> cum_inc;
  vector<Type> surv;
  
  // Constructor (just a dummy! Would like to change to be model indicator)
  report_values(SEXP x){ 
    u_fixed_mmc = asVector<Type>(getListElement(x,"u_fixed_mmc"));
  }
};
