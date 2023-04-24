/// @file threemc.h
// #include <TMB.hpp>

// Data Struct //
#ifndef THREEMC_DATA_DEF
#define THREEMC_DATA_DEF

template<class Type>
struct Threemc_data {

  // indicators
  int is_type;         // Model with type (MMC/TMC) split, or without?
  int rw_order;        // Model with AR1 or RW temporal prior?
  int paed_age_cutoff; // Model with paedaitric age cutoff for medical circumcisions?
  int inc_time_tmc;    // Model with time TMC effect?

  // Survival analysis matrices
  density::SparseMatrix<Type> A_mmc; // Matrix selecting instantaneous hazard for medically circumcised pop
  density::SparseMatrix<Type> A_tmc; // Matrix selecting instantaneous hazard for traditionally circumcised pop
  density::SparseMatrix<Type> A_mc; // Matrix selecting instantaneous hazard for unknown circumcised pop
  density::SparseMatrix<Type> B; // Matrix selecting relevant cumulative hazard entry for observed and right censored pop
  density::SparseMatrix<Type> C; // Matrix selecting relevant cumulative hazard entry for interval censored pop

  // Integeration matrices
  density::SparseMatrix<Type> IntMat1;
  density::SparseMatrix<Type> IntMat2;

  // Design Matrices
  density::SparseMatrix<Type> X_fixed_mmc; // Design matrix for the fixed effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_time_mmc; // Design matrix for the temporal random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_age_mmc; // Design matrix for the stratification random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_space_mmc; // Design matrix for the stratification random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_agetime_mmc; // Design matrix for the interaction random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_spacetime_mmc; // Design matrix for the interaction random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_agespace_mmc; // Design matrix for the interaction random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_fixed_tmc; // Design matrix for the fixed effects in the traditional circumcision hazard rate
  density::SparseMatrix<Type> X_age_tmc; // Design matrix for the stratification random effects in the traditional circumcision hazard rate
  density::SparseMatrix<Type> X_space_tmc; // Design matrix for the stratification random effects in the medical circumcision hazard rate
  density::SparseMatrix<Type> X_agespace_tmc; // Design matrix for the interaction random effects in the medical circumcision hazard rate

  // For model with paediatric age cutoff
  density::SparseMatrix<Type> X_fixed_mmc_paed; // Design matrix for the fixed effects in the medical circumcision hazard rate (under specified age)
  density::SparseMatrix<Type> X_age_mmc_paed; 
  density::SparseMatrix<Type> X_space_mmc_paed; 
  density::SparseMatrix<Type> X_agespace_mmc_paed;

  // For model with time TMC effect
  density::SparseMatrix<Type> X_time_tmc; // Design matrix for the temporal random effects in the traditional circumcision hazard rate
  density::SparseMatrix<Type> X_agetime_tmc; // Design matrix for the interaction random effects in the traditional circumcision hazard rate
  density::SparseMatrix<Type> X_spacetime_tmc; // Design matrix for the interaction random effects in the traditional circumcision hazard rate

  // for model with no type
  density::SparseMatrix<Type> X_fixed;    // Design matrix for the fixed effects
  density::SparseMatrix<Type> X_time;     // Design matrix for the temporal random effects
  density::SparseMatrix<Type> X_age;      // Design matrix for the stratification random effects
  density::SparseMatrix<Type> X_space;    // Design matrix for the stratification random effects
  density::SparseMatrix<Type> X_agetime;  // Design matrix for the age-time interaction random effects
  density::SparseMatrix<Type> X_agespace; // Design matrix for the age-space interaction random effects
  density::SparseMatrix<Type> X_spacetime; // Design matrix for the age-space interaction random effects
 
  // Aggr mat for # circs (if RW, precision mat for spatial REs)
  density::SparseMatrix<Type> Q_space; 
  density::SparseMatrix<Type> Q_time;  // Precision matrix for spatial process (only for RW model)

  // Constructor
  Threemc_data(SEXP x);
};

#endif


/*******************************************************************/
/*   Class to calculate negative log likelihood in threemc model   */
/*******************************************************************/

#ifndef THREEMC_DEF
#define THREEMC_DEF

using namespace density;

// Model with type (MMC/TMC) split, no paed age cutoff or time TMC effect
template <class Type>
class Threemc {

  protected:

    // negative log likelihood
    Type nll; 
  
    // report values (hazard rates, incidence and cumulative incidence)
    vector<Type> haz_mmc;     // Medical hazard rate
    vector<Type> haz_tmc;     // Traditional hazard rate
    vector<Type> haz;         // Total hazard rate
    vector<Type> inc_tmc;     // Traditional circumcision incidence rate
    vector<Type> inc_mmc;     // Medical circumcision incidence rate
    vector<Type> inc;         // Total circumcision incidence rate
    vector<Type> cum_inc_tmc; // Traditional circumcision cumulative incidence rate
    vector<Type> cum_inc_mmc; // Medical circumcision cumulative incidence rate
    vector<Type> cum_inc;     // Total circumcision cumulative incidence rate
    vector<Type> surv;        // Survival probabilities
    vector<Type> surv_lag;    // Lagged survival probabilities
    vector<Type> leftcens;    // left censored incidence rate

  public:

    // Default Constructor
    Threemc();

    // Default virtual Destructor
    virtual ~Threemc();

	  // Fixed effects for circumcision rate
    void fix_eff_p(vector<Type> u_fixed);

    // Prior on temporal random effects
    void rand_eff_time_p(vector<Type> u_time,
                         Type logsigma_time,
                         Type sigma_time,
                         Type logitrho_time1,
                         Type rho_time1);

    // Prior on the age random effects
    void rand_eff_age_p(vector<Type> u_age,
                        Type logsigma_age,
                        Type sigma_age,
                        Type logitrho_age1,
                        Type rho_age1);

    // Prior on the spatial random effects
    void rand_eff_space_p(density::SparseMatrix<Type> Q_space,
                          vector<Type> u_space,
                          Type logsigma_space,
                          Type sigma_space);

    // Sum to zero constraints
    void sum_to_zero(array<Type> u_interact);

    // Prior on the interaction random effects for either MMC or MC (for model w/ no type)
    void rand_eff_interact_p(density::SparseMatrix<Type> Q_space,
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
                             Type rho_time3);

    // Overload prior on the interaction random effects for just TMC
    void rand_eff_interact_p(density::SparseMatrix<Type> Q_space,
                             array<Type> u_agespace_tmc,
                             Type logsigma_agespace_tmc,
                             Type sigma_agespace_tmc,
                             Type logitrho_tmc_age2,
                             Type rho_tmc_age2);

    // Function to calculate report values 
    // TODO: This will change depending on whether type information is included
    // Need to just overload this function
    // TODO: Can definitely design this function better to avoid repitition
    // Can pass reference to haz to do this (should work!)
    // For MMC: 
    void calc_haz(vector<Type> &hazard,
                  density::SparseMatrix<Type> X_fixed, 
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
                  Type sigma_spacetime,
                  int scale, // should we scale hazard on [0, 1]?
                  int init); // should we initialise or add to hazard?

    // For TMC: 
  void calc_haz(vector<Type> &hazard,
                density::SparseMatrix<Type> X_fixed,
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
                Type sigma_agespace,
                int scale,
                int init);

    // final calculation of total report vals (e.g. haz = haz_mmc + haz_tmc)
    void calc_haz();

    // function to calculate survival probabilities
    void calc_surv(density::SparseMatrix<Type> IntMat1,
                   density::SparseMatrix<Type> IntMat2);

    // Function to calculate incidence & cumulative incidence
    // TODO: This code is reused, find a way to change e.g inc_mmc while referring to inc here
    // Actually works well here with if statement, but may not scale well with more models
    void calc_inc(density::SparseMatrix<Type> IntMat1, int is_type);
    // void calc_inc(density::SparseMatrix<Type> IntMat1);
    // void calc_inc(density::SparseMatrix<Type> IntMat1, int is_type);

    // Function to calculate likelihood
    // TODO: Can make an enum (or something?) pointer to iterate over for this
    // (will be different for each model)
    void likelihood(density::SparseMatrix<Type> Mat,
                    vector<Type> report_val);

    void calc_nll(struct Threemc_data<Type> threemc_data, objective_function<Type>* obj);

    // Getter for nll
    Type get_nll() {
      return nll;
    };
};
#endif

// Model with type split and random walk (not autoregressive) temporal prior

#ifndef THREEMC_RW_DEF
#define THREEMC_RW_DEF

// Model with no type split, no paed age cutoff or time TMC effect
template<class Type>
class Threemc_rw : virtual public Threemc<Type> {

  protected:

    using Threemc<Type>::nll;
    using Threemc<Type>::haz_mmc;
    using Threemc<Type>::haz_tmc;
    using Threemc<Type>::haz;
    using Threemc<Type>::inc_mmc;
    using Threemc<Type>::inc_tmc;
    using Threemc<Type>::inc;
    using Threemc<Type>::cum_inc_mmc;
    using Threemc<Type>::cum_inc_tmc;
    using Threemc<Type>::cum_inc;
    using Threemc<Type>::surv;
    using Threemc<Type>::surv_lag;
    using Threemc<Type>::leftcens;
 
  public:

    // Default Constructor
    Threemc_rw();

    // Default virtual Destructor
    virtual ~Threemc_rw();

    // Base functions
    // TODO: Remove functions we overload here!
    using Threemc<Type>::fix_eff_p;
    using Threemc<Type>::rand_eff_age_p;
    using Threemc<Type>::rand_eff_space_p;
    using Threemc<Type>::sum_to_zero;
    using Threemc<Type>::rand_eff_interact_p;
    using Threemc<Type>::calc_haz;
    using Threemc<Type>::calc_surv;
    using Threemc<Type>::calc_inc;
    using Threemc<Type>::likelihood;
    using Threemc<Type>::get_nll;

    void rand_eff_time_p(density::SparseMatrix<Type> Q_time,
                         vector<Type> u_time,
                         Type logsigma_time,
                         Type sigma_time);

    void rand_eff_interact_p(density::SparseMatrix<Type> Q_space,
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
                             Type rho_age3);

    void calc_nll(struct Threemc_data<Type> threemc_data,
                  objective_function<Type>* obj);
};

#endif

#ifndef THREEMC_PAED_DEF
#define THREEMC_PAED_DEF

// Model with paediatric age cutoff for MMC
template<class Type>
class Threemc_paed : virtual public Threemc<Type> {

  protected:

    using Threemc<Type>::nll;
    using Threemc<Type>::haz_mmc;
    using Threemc<Type>::haz_tmc;
    using Threemc<Type>::haz;
    using Threemc<Type>::inc_mmc;
    using Threemc<Type>::inc_tmc;
    using Threemc<Type>::inc;
    using Threemc<Type>::cum_inc_mmc;
    using Threemc<Type>::cum_inc_tmc;
    using Threemc<Type>::cum_inc;
    using Threemc<Type>::surv;
    using Threemc<Type>::surv_lag;
    using Threemc<Type>::leftcens;
 
  public:

    // Default Constructor
    Threemc_paed();

    // Default virtual Destructor
    virtual ~Threemc_paed();

    // Base functions
    using Threemc<Type>::fix_eff_p;
    using Threemc<Type>::rand_eff_time_p;
    using Threemc<Type>::rand_eff_age_p;
    using Threemc<Type>::rand_eff_space_p;
    using Threemc<Type>::sum_to_zero;
    using Threemc<Type>::rand_eff_interact_p;
    using Threemc<Type>::calc_haz; // only using TMC version of this function
    using Threemc<Type>::calc_surv;
    using Threemc<Type>::calc_inc;
    using Threemc<Type>::likelihood;
    using Threemc<Type>::get_nll;

    void calc_nll(struct Threemc_data<Type> threemc_data,
                  objective_function<Type>* obj);
};

#endif


#ifndef THREEMC_TIME_TMC_DEF
#define THREEMC_TIME_TMC_DEF

// Model with time effect for TMC
template<class Type>
class Threemc_time_tmc : virtual public Threemc<Type> {

  protected:

    using Threemc<Type>::nll;
    using Threemc<Type>::haz_mmc;
    using Threemc<Type>::haz_tmc;
    using Threemc<Type>::haz;
    using Threemc<Type>::inc_mmc;
    using Threemc<Type>::inc_tmc;
    using Threemc<Type>::inc;
    using Threemc<Type>::cum_inc_mmc;
    using Threemc<Type>::cum_inc_tmc;
    using Threemc<Type>::cum_inc;
    using Threemc<Type>::surv;
    using Threemc<Type>::surv_lag;
    using Threemc<Type>::leftcens;
 
  public:

    // Default Constructor
    Threemc_time_tmc();

    // Default virtual Destructor
    virtual ~Threemc_time_tmc();

    // Base functions
    using Threemc<Type>::fix_eff_p;
    using Threemc<Type>::rand_eff_age_p;
    using Threemc<Type>::rand_eff_time_p; // run for TMC as for MMC
    using Threemc<Type>::rand_eff_space_p;
    using Threemc<Type>::sum_to_zero;
    using Threemc<Type>::rand_eff_interact_p;
    using Threemc<Type>::calc_haz;
    using Threemc<Type>::calc_surv;
    using Threemc<Type>::calc_inc;
    using Threemc<Type>::likelihood;
    using Threemc<Type>::get_nll;

    // Need calc_haz where there is a time effect but no time interactions
    void calc_haz(vector<Type> &hazard,
                  density::SparseMatrix<Type> X_fixed, 
                  density::SparseMatrix<Type> X_age, 
                  density::SparseMatrix<Type> X_time, 
                  density::SparseMatrix<Type> X_space,
                  density::SparseMatrix<Type> X_agespace,
                  // parameters
                  vector<Type> u_fixed,
                  vector<Type> u_age,
                  vector<Type> u_time,
                  vector<Type> u_space,
                  array<Type> u_agespace,
                  Type sigma_age,
                  Type sigma_time,
                  Type sigma_space,
                  Type sigma_agespace,
                  int scale,
                  int init);

    void calc_nll(struct Threemc_data<Type> threemc_data,
                  objective_function<Type>* obj);
};

#endif


// Model with paediatric age cutoff for MMC and a RW temporal prior
#ifndef THREEMC_PAED_RW_DEF
#define THREEMC_PAED_RW_DEF

template<class Type>
class Threemc_paed_rw : virtual public Threemc_rw<Type>, virtual public Threemc_paed<Type> {

  protected:

    using Threemc<Type>::nll;
    using Threemc<Type>::haz_mmc;
    using Threemc<Type>::haz_tmc;
    using Threemc<Type>::haz;
    using Threemc<Type>::inc_mmc;
    using Threemc<Type>::inc_tmc;
    using Threemc<Type>::inc;
    using Threemc<Type>::cum_inc_mmc;
    using Threemc<Type>::cum_inc_tmc;
    using Threemc<Type>::cum_inc;
    using Threemc<Type>::surv;
    using Threemc<Type>::surv_lag;
    using Threemc<Type>::leftcens;
 
  public:

    // Default Constructor
    Threemc_paed_rw();

    // Default virtual Destructor
    virtual ~Threemc_paed_rw();

    // Base functions
    using Threemc<Type>::fix_eff_p;
    using Threemc_rw<Type>::rand_eff_time_p;
    using Threemc<Type>::rand_eff_age_p;
    using Threemc<Type>::rand_eff_space_p;
    using Threemc<Type>::sum_to_zero;
    using Threemc_rw<Type>::rand_eff_interact_p;
    using Threemc<Type>::calc_haz; 
    using Threemc_paed<Type>::calc_haz; 
    using Threemc<Type>::calc_surv;
    using Threemc<Type>::calc_inc;
    using Threemc<Type>::likelihood;
    using Threemc<Type>::get_nll;

    void calc_nll(struct Threemc_data<Type> threemc_data,
                  objective_function<Type>* obj);
};

#endif


// Model with time effect for TMC and paediatric age cutoff for MMC 
#ifndef THREEMC_PAED_TIME_TMC_DEF
#define THREEMC_PAED_TIME_TMC_DEF

template<class Type>
class Threemc_paed_time_tmc : virtual public Threemc_paed<Type>,
                              virtual public Threemc_time_tmc<Type> {

  protected:

    using Threemc<Type>::nll;
    using Threemc<Type>::haz_mmc;
    using Threemc<Type>::haz_tmc;
    using Threemc<Type>::haz;
    using Threemc<Type>::inc_mmc;
    using Threemc<Type>::inc_tmc;
    using Threemc<Type>::inc;
    using Threemc<Type>::cum_inc_mmc;
    using Threemc<Type>::cum_inc_tmc;
    using Threemc<Type>::cum_inc;
    using Threemc<Type>::surv;
    using Threemc<Type>::surv_lag;
    using Threemc<Type>::leftcens;
 
  public:

    // Default Constructor
    Threemc_paed_time_tmc();

    // Default virtual Destructor
    virtual ~Threemc_paed_time_tmc();

    // Base functions
    using Threemc<Type>::fix_eff_p;
    using Threemc<Type>::rand_eff_age_p;
    using Threemc<Type>::rand_eff_time_p; // run for TMC as for MMC
    using Threemc<Type>::rand_eff_space_p;
    using Threemc<Type>::sum_to_zero;
    using Threemc<Type>::rand_eff_interact_p;
    using Threemc<Type>::calc_haz;
    using Threemc_time_tmc<Type>::calc_haz;
    using Threemc<Type>::calc_surv;
    using Threemc<Type>::calc_inc;
    using Threemc<Type>::likelihood;
    using Threemc<Type>::get_nll;

    void calc_nll(struct Threemc_data<Type> threemc_data,
                  objective_function<Type>* obj);
};

#endif


// Model with no type split, no paed age cutoff or time TMC effect
#ifndef THREEMC_NT_DEF
#define THREEMC_NT_DEF

template<class Type>
class Threemc_nt : virtual public Threemc<Type> {

  protected:

    using Threemc<Type>::nll;
    using Threemc<Type>::haz;
    using Threemc<Type>::inc;
    using Threemc<Type>::cum_inc;
    using Threemc<Type>::surv;
    using Threemc<Type>::surv_lag;
    using Threemc<Type>::leftcens;
 
  public:

    // Default Constructor
    Threemc_nt();

    // Default virtual Destructor
    virtual ~Threemc_nt();

    // Base functions
    using Threemc<Type>::fix_eff_p;
    using Threemc<Type>::rand_eff_time_p;
    using Threemc<Type>::rand_eff_age_p;
    using Threemc<Type>::rand_eff_space_p;
    using Threemc<Type>::sum_to_zero;
    using Threemc<Type>::rand_eff_interact_p;
    using Threemc<Type>::calc_haz;
    using Threemc<Type>::calc_surv;
    using Threemc<Type>::calc_inc;
    using Threemc<Type>::likelihood;
    using Threemc<Type>::get_nll;
 
    void calc_nll(struct Threemc_data<Type> threemc_data,
                  objective_function<Type>* obj);
};

#endif

#ifndef THREEMC_NT_RW_DEF
#define THREEMC_NT_RW_DEF

// Model with no type split which uses a RW temporal prior
template<class Type>
class Threemc_nt_rw : virtual public Threemc_nt<Type>, virtual public Threemc_rw<Type> {

  protected:

    using Threemc<Type>::nll;
    using Threemc<Type>::haz;
    using Threemc<Type>::inc;
    using Threemc<Type>::cum_inc;
    using Threemc<Type>::surv;
    using Threemc<Type>::surv_lag;
    using Threemc<Type>::leftcens;
 
  public:

    // Default Constructor
    Threemc_nt_rw();

    // Default virtual Destructor
    virtual ~Threemc_nt_rw();

    // Base functions
    using Threemc<Type>::fix_eff_p;
    using Threemc<Type>::rand_eff_time_p;
    using Threemc_rw<Type>::rand_eff_time_p;
    using Threemc<Type>::rand_eff_age_p;
    using Threemc<Type>::rand_eff_space_p;
    using Threemc<Type>::sum_to_zero;
    using Threemc_rw<Type>::rand_eff_interact_p;
    using Threemc_nt<Type>::calc_haz;
    using Threemc<Type>::calc_surv;
    using Threemc<Type>::calc_inc;
    using Threemc<Type>::likelihood;
    using Threemc<Type>::get_nll;

    void calc_nll(struct Threemc_data<Type> threemc_data,
                  objective_function<Type>* obj);
};

#endif

