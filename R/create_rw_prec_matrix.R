
#' @title Create Precision Matrix for RWp Process
#'
#' @description Create the precision matrix for a RWp process.
#'
#' @param dim Dimension of the precision matrix.
#' @param order Order of the random walk, Default: 1
#' @param offset.diag Option to offset diagonal by 1E-6, Default: TRUE
#'
#' @seealso
#'  \code{\link[methods]{as}}
#
#' @return RW precision matrix
#' @export
create_rw_prec_matrix <- function(dim,
                                  order = 1,
                                  offset.diag = TRUE) {
  ## Creating structure matrix
  Q <- diff(diag(dim), differences = order)
  Q <- t(Q) %*% Q
  ## Adding offset to diagonal if required
  if (offset.diag) {
    diag(Q) <- diag(Q) + 1E-6
  }
  ## Converting to sparse matrix
  Q <- methods::as(Q, "sparseMatrix")
  ## Returning matrix
  return(Q)
}
