#' @title Compile TMB model without C++ Warnings
#' @description Wrapper around \code{\link[TMB]{compile}} which optionally 
#' redirects C++ compiler output to a logfile, avoiding your R console, since 
#' it can be quite verbose. Note that \code{logfile} argument default is not
#' valid for Windows.
#' @param file C++ file.
#' @param logfile File to redirect output from C++ compiler. To print output 
#' to R console, specify \code{logfile = NULL}, Default: '/tmp/tmb_logfile.log'
#' @param ... Additional parameters for \code{\link[TMB]{compile}}.
#' 
#' @seealso 
#'  \code{\link[TMB]{compile}}
#' @rdname compile_tmb
#' @export
compile_tmb <- function(file, 
                        logfile = "/tmp/tmb_logfile.log",
                        ...) {

  if (!is.null(logfile)) {

    if (!file.exists(logfile)) file.create(logfile)
    logfile_redirect <- paste0("&> ", logfile)
    
    tryCatch({
      invisible(TMB::compile(file, logfile_redirect, ...))
      message(
        paste("any output from TMB::compile has been redirected to \n",
              logfile, 
              "\n Please specify 'logfile = NULL' to print to your R console")
      )
    }, error = function(e) {
      stop("TMB::compile has produced an error.\n
           Please specify 'logfile = NULL' to return function error messages")
    }) 
  } else {
    tmb::compile(file, ...)
  }
}
