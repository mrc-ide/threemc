.onLoad <- function(libname, pkgname) {
  message("Setting spdep::set.ZeroPolicyOption(TRUE)")
  spdep::set.ZeroPolicyOption(TRUE)
}

