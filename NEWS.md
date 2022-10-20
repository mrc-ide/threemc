# threemc 0.1.26

* Fixed bug in `prepare_survey_data` whereby surveys not at maximum area 
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

