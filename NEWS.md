# threemc 0.1.1

* Added a `NEWS.md` file to track changes to the package.
* Change title: (Matt's) Multi-Level Model of Male Circumcision in Sub-Saharan Africa
* Remove wildcard matching (`dplyr::any_of()`, `dplyr::contains()`) from survey data 
  set construction.
* Don't reference the `space` index in the survey data set; only add when merging 
  survey data to the model frame.
  
# threemc 0.1.2 

* Conditionals added to produce design matrices when modelling at the country 
level, where the adjacency matrix will just be a 1x1 matrix with entry 0. 
