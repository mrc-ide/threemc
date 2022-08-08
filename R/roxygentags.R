#' @title Export TMB models
#'
#' This is a dummy function who's purpose is to hold the useDynLib roxygen tag.
#' This tag will populate the namespace with compiled c++ functions upon
#' package install.
#'
#' @useDynLib Surv_SpaceAgeTime
#' @useDynLib Surv_SpaceAgeTime_ByType_withUnknownType
#' @useDynLib Surv_SpaceAgeTime_ByType_withUnknownType_No_AgeTime_Interaction
#' @noRd
#' @keywords internal
.dummy <- function() {
  return(NULL)
}
