# threemc 0.1.38

- Fix bug in `threemc_ppc`. Weighted mean is calculated using `circ_status` and 
`indweight`. However, this ignores `circ_type`. For `type %in% c("MMC", "TMC")`,
need to have `circ_status == 0` for all `circ_type == "Missing"`. This does not 
apply to `type == "MC"`, as we do not need to know `circ_type` to still count it 
amongst `MC` circumcisions.

# threec 0.1.37

- In `create_shell_dataset`, replace missing populations for less granular areas
with their respective child areas' populations. 

# threemc 0.1.36

- Add optional penalised time spline, by specifying a non-null value for the 
`k_dt_time` argument of `threemc_prepare_model_data`.

# threemc 0.1.35

- Adds models 
("Surv_SpaceAgeTime_ByType_withUnknownType_RW_MMC2" and 
"Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW_MMC2") which use 
a random walk temporal prior for MMC, but an AR 1 temporal prior for TMC.


# threemc 0.1.34 

- Version 0.1.27 introduced filling in missing population data in 
`create_shell_dataset`; this update decouples this behaviour from this function 
into a separate internal function `fill_downup_populations`, and this 
functionality has been added to `threemc_aggregate` as well. 

# threemc 0.1.33

- `threemc_ppc` rewritten in `data.table` for significant speed and memory 
efficicency increase.

# threemc 0.1.32 

- Add initial introductory vignette.

# threemc 0.1.31

* Add TMB model for every iteration of: 
  - Including a time effect for TMC or not, 
  - Including a constant peadiatric MMC rate or not, and 
  - Using a random walk, rather than AR 1 prior
Urgent need to functionalise all of this to greatly decrease code complexity 
and package compilation time. 

# threemc 0.1.30

* Have verbose output on function progress from `threemc_fit_model` as default.

# threemc 0.1.29

* Can have `threemc_fit_model` choose `mod` itself based on parameters in 
either `parameters` (for new fits) or `fit$par` (for re-sampling from fits). 
Like 0.1.28, this update abstracts model specification from the end user, which 
involves generally involves setting a global variable `mod` to a quite 
"esoteric" model name, such as "Surv_SpaceAgeTime_ByType_withUnknownType_RW2". 
Instead, the user can simply specify the much more intuitive parameters 
`rw_order` and `include_tmc` to `threemc_prepare_model_data`, simplifying their 
experience.


# threemc 0.1.28

* new function `threemc_initial_pars` to abstract initial hyperparameter 
specification from the end user. In scripts, this section is quite long and 
ugly and in general can often go wrong due to parameter order etc, so it is 
best to abstract this functionality. Defaults can still be overridden (e.g. 
when we want to fit a model with "mapped"/fixed hyperparameters) by specifying 
the `custom_init` argument using a named list of parameter values.


# threemc 0.1.27 

* In `prepare_survey_data`, fill `NA` populations for earlier years than we 
have data for with the earliest known value for each `age` and `area_id`. 

# threemc 0.1.26

* Fix bug in `prepare_survey_data` whereby surveys not at maximum area 
level have their `area_id` columns coerced to NA, when reassiging 
survey `area_level` to the specified argument `area_level`.


# threemc 0.1.25

* Add models which include a random effect for time for traditional 
male circumcision. 
* Add `inc_time_tmc` argument to `threemc_prepare_model_data` to 
produce `X_time_tmc` so we can use these non-constant TMC models. 

# threemc 0.1.24

* Add function (`threemc_oos_pcc`) to perform posterior predictive checks 
with OOS survey estimates, for model validation and comparison. 

# threemc 0.1.23

* Add function (`survey_points_dmmpt2_convert_convention`) which can be used 
to change the age group, circumcision type and column naming conventions 
used in DMMPT2 and empirical survey circumcision estimates to match those 
used in threemc aggregations. 

# threemc 0.1.22

* Add `rw_order` argument to `threemc_prepare_model_data`, which allows 
one to specify a Random Walk temporal process for our temporal prior. 
Leaving `rw_order = NULL` uses the default AR 1 temporal prior.

# threemc 0.1.21

* Replace loop in `aggregate_sample_age_group` with method which uses 
`left_join`, ~4x faster for SWZ, likely scales much better for larger countries. 

# threemc 0.1.2 

* Conditionals added to produce design matrices when modelling at the country 
level, where the adjacency matrix will just be a 1x1 matrix with entry 0. 

# threemc 0.1.1

* Added a `NEWS.md` file to track changes to the package.
* Change title: (Matt's) Multi-Level Model of Male Circumcision in Sub-Saharan Africa
* Remove wildcard matching (`dplyr::any_of()`, `dplyr::contains()`) from survey data 
  set construction.
* Don't reference the `space` index in the survey data set; only add when merging 
  survey data to the model frame.

