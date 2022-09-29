#define TMB_LIB_INIT R_init_threemc
#include <TMB.hpp>

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
/************************************************************************/
/* Objective function to specify model and to optimize model parameters */
/************************************************************************/
template<class Type>
Type objective_function<Type>::operator() ()
{
	using namespace density;

    ////////////////////////
    /// Data definitions ///
    ////////////////////////
	// Survival analysis matrices
    DATA_SPARSE_MATRIX(A1); // Matrix selecting instantaneous hazard for circumcised pop
    DATA_VECTOR(A2); // Weighting for relevant instantaneous hazard for circumcised pop in the likelihood
    DATA_SPARSE_MATRIX(B1); // Matrix selecting relevant cumulative hazard entry for observed and right censored pop
    DATA_VECTOR(B2); // Weighting for relevant cumulative hazard entry for observed and right censored pop in the likelihood
    DATA_SPARSE_MATRIX(C1); // Matrix selecting relevant cumulative hazard entry for interval censored pop
    DATA_VECTOR(C2); // // Weighting for relevant cumulative hazard entry for interval censored pop in the likelihood
    DATA_SPARSE_MATRIX(IntMat1); // Integration matrix for cumulative hazard
    DATA_SPARSE_MATRIX(IntMat2); // Integration matrix for lagged cumulative hazard

	// Design matrices
    DATA_SPARSE_MATRIX(X_fixed);    // Design matrix for the fixed effects
    DATA_SPARSE_MATRIX(X_time);     // Design matrix for the temporal random effects
    DATA_SPARSE_MATRIX(X_age);      // Design matrix for the stratification random effects
    DATA_SPARSE_MATRIX(X_space);    // Design matrix for the stratification random effects
    DATA_SPARSE_MATRIX(X_agetime);  // Design matrix for the age-time interaction random effects
    DATA_SPARSE_MATRIX(X_agespace); // Design matrix for the age-space interaction random effects
    DATA_SPARSE_MATRIX(X_spacetime); // Design matrix for the age-space interaction random effects

	// Adjacency matrices
    DATA_SPARSE_MATRIX(Q_space); // Precision matrix for spatial random effects
	DATA_SPARSE_MATRIX(Q_time); // Precision matrix for spatial process

    //////////////////
    /// Parameters ///
    //////////////////
	// Fixed effects
    PARAMETER_VECTOR(u_fixed);

	// Age random effects
    PARAMETER_VECTOR(u_age);

	// Temporal random effects
    PARAMETER_VECTOR(u_time);

	// Spatial random effects
    PARAMETER_VECTOR(u_space);

	// Interaction random effects
    PARAMETER_ARRAY(u_agetime);
    PARAMETER_ARRAY(u_agespace);
    PARAMETER_ARRAY(u_spacetime);

	// Standard deviations
    PARAMETER(logsigma_age);       Type sigma_age       = exp(logsigma_age);
    PARAMETER(logsigma_time);      Type sigma_time      = exp(logsigma_time);
    PARAMETER(logsigma_space);     Type sigma_space     = exp(logsigma_space);
    PARAMETER(logsigma_agetime);   Type sigma_agetime   = exp(logsigma_agetime);
    PARAMETER(logsigma_agespace);  Type sigma_agespace  = exp(logsigma_agespace);
    PARAMETER(logsigma_spacetime); Type sigma_spacetime = exp(logsigma_spacetime);

	// Autocorrelation parameters
    PARAMETER(logitrho_age1);  Type rho_age1   = geninvlogit(logitrho_age1,  Type(-1.0), Type(1.0));
    PARAMETER(logitrho_age2);  Type rho_age2   = geninvlogit(logitrho_age2,  Type(-1.0), Type(1.0));
    PARAMETER(logitrho_age3);  Type rho_age3   = geninvlogit(logitrho_age3,  Type(-1.0), Type(1.0));

	//////////////////////////////////
	/// Prior on the fixed effects ///
	//////////////////////////////////
	// Negative log likelihood definition
    Type nll = Type(0);

	// Fixed effects for the circumcision rate
    nll -= sum(dnorm(u_fixed,  Type(0), Type(5), TRUE));

	////////////////////////////////////////////
	/// Prior on the temporal random effects ///
	////////////////////////////////////////////
    // AR1 Process 
    nll += GMRF(Q_time)(u_time);
  
    // Sum to zero constraint 
    nll -= dnorm(u_time.sum(), Type(0), Type(0.001) * u_time.size(), TRUE);
  
    // Prior on the standard deviation for the temporal random effects
    nll -= dexp(sigma_time, Type(1), TRUE) + logsigma_time;

	///////////////////////////////////////
	/// Prior on the age random effects ///
	///////////////////////////////////////
	// AR1 Process
	nll += AR1(rho_age1)(u_age);

	// Sum to zero constraint
	nll -= dnorm(u_age.sum(), Type(0), Type(0.001) * u_age.size(), TRUE);

    // Prior on the standard deviation for the age random effects
	nll -= dexp(sigma_age, Type(1), TRUE) + logsigma_age;

	// Prior on the logit autocorrelation parameters
	nll -= dnorm(logitrho_age1, Type(3), Type(1), TRUE);

	//////////////////////////////////////////////////
	/// Prior on the stratification random effects ///
	//////////////////////////////////////////////////
	// Gaussian markov random field with prespecified precision matrix
	nll += GMRF(Q_space)(u_space);

	// Sum to zero constraint
	nll -= dnorm(u_space.sum(), Type(0), Type(0.001) * u_space.size(), TRUE);

    // Prior on the standard deviation for the stratification random effects
	nll -= dexp(sigma_space, Type(1), TRUE) + logsigma_space;

	///////////////////////////////////////////////
	/// Prior on the interaction random effects ///
	///////////////////////////////////////////////
	// Interactions: space-time (GMRF x AR1), age-time (AR1 x AR1) and age-space (AR1 x GMRF)
	nll += SEPARABLE(GMRF(Q_time), AR1(rho_age2))(u_agetime);
	nll += SEPARABLE(GMRF(Q_space), AR1(rho_age3))(u_agespace);
	nll += SEPARABLE(GMRF(Q_space), GMRF(Q_time))(u_spacetime);

    // Sum-to-zero constraints
    for (int i = 0; i < u_agetime.cols(); i++) {
      nll -= dnorm(u_agetime.col(i).sum(), Type(0), Type(0.001) * u_agetime.col(i).size(), true);
  	}
	for (int i = 0; i < u_agespace.cols(); i++) {
      nll -= dnorm(u_agespace.col(i).sum(), Type(0), Type(0.001) * u_agespace.col(i).size(), true);
  	}
	for (int i = 0; i < u_spacetime.cols(); i++) {
      nll -= dnorm(u_spacetime.col(i).sum(), Type(0), Type(0.001) * u_spacetime.col(i).size(), true);
  	}

	// Vectorising the interaction
	vector<Type> u_agetime_v(u_agetime);
	vector<Type> u_agespace_v(u_agespace);
	vector<Type> u_spacetime_v(u_spacetime);

	// Prior on the standard deviations
	nll -= dexp(sigma_agetime, Type(1), TRUE) + logsigma_agetime;
	nll -= dexp(sigma_agespace, Type(1), TRUE) + logsigma_agespace;
	nll -= dexp(sigma_spacetime, Type(1), TRUE) + logsigma_spacetime;

	// Prior on the logit autocorrelation parameters
	nll -= dnorm(logitrho_age2, Type(3), Type(1), TRUE);
	nll -= dnorm(logitrho_age3, Type(3), Type(1), TRUE);

    //////////////////////////////
    /// Estimating hazard rate ///
    //////////////////////////////
    // Hazard Rate
    vector<Type> haz = X_fixed * u_fixed +
		X_age * u_age * sigma_age +
			X_space * u_space * sigma_space +
				X_time * u_time * sigma_time +
					X_agetime * u_agetime_v * sigma_agetime +
						X_agespace * u_agespace_v * sigma_agespace +
							X_spacetime * u_spacetime_v * sigma_spacetime;

	// Rates on [0,1] scale
    haz = invlogit_vec(haz);

	// Survival probabilities
	vector<Type> logprob  = log(Type(1.0) - haz);
	vector<Type> surv     = exp(IntMat1 * logprob);
	vector<Type> surv_lag = exp(IntMat2 * logprob);

	// Incidence
	vector<Type> inc = haz * surv_lag;

	// Cumulative incidence
	vector<Type> cum_inc = IntMat1 * inc;

    //////////////////
    /// Likelihood ///
    //////////////////
	// Getting likelihood for those circumcised
	nll -= (A2 * log(A1 * inc)).sum();

	// Getting likelihood for those right censored
	nll -= (B2 * log(B1 * surv)).sum();

	// Getting likelihood for those left censored
	nll -= (C2 * log(C1 * cum_inc)).sum();

    ///////////////////////////
    /// Reporting variables ///
    ///////////////////////////
    REPORT(haz);     // Hazard rate
    REPORT(surv);    // Survival probabilities
    REPORT(inc);     // Circumcision incidence rate
    REPORT(cum_inc); // Circumcision cumulative incidence rate

    /////////////////////////////////////////
    /// Returning negative log likelihood ///
    /////////////////////////////////////////
    return nll;
}
