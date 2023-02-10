## check package with goodpractise ##

stopifnot(basename(getwd()) == "threemc") # fail if not running for threemc

# assume you don't want to build vignettes; very slow!
if (!exists(check_vignettes) && !rlang::is_bool(check_vignettes)) {
  check_vignettes <- FALSE
}

# run devtools::check() instead of rcmdcheck()
checks <- goodpractice::all_checks()[
  !grepl("rcmdcheck", goodpractice::all_checks())
]
goodpractice::goodpractice(
  checks = checks,
  extra_checks = devtools::check(
    vignettes = check_vignettes, # set to TRUE if changing vignettes
    build_args = "--resave-data" # build with compressed datasets
  ),
  quiet = FALSE
)
