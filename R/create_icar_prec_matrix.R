#' @title Create Precision Matrix for ICAR Process
#'
#' @description Create the precision matrix for an ICAR process.
#'
#' @param sf_obj Shapefiles needed for adjacency.
#' @param row.names Unique IDs for the areas.
#' @return ICAR precision matrix.
#'
#' @seealso
#'  \code{\link[spdep]{poly2nb}}
#'  \code{\link[spdep]{nb2mat}}
#'  \code{\link[naomi]{scale_gmrf_precision}}
#' @rdname create_icar_prec_matrix
#' @export
create_icar_prec_matrix <- function(sf_obj = NULL,
                                    row.names = NULL) {
  ## Creating neighbourhood structure
  Q_space <- spdep::poly2nb(sf_obj, row.names = sf_obj[, row.names])
  ## Converting to adjacency matrix
  Q_space <- spdep::nb2mat(Q_space, style = "B", zero.policy = TRUE)
  ## Converting to sparse matrix
  Q_space <- as(Q_space, "sparseMatrix")
  ## Creating precision matrix from adjacency
  Q_space <- naomi::scale_gmrf_precision(
      Q   = diag(rowSums(as.matrix(Q_space))) - 0.99 * Q_space,
      A   = matrix(1, 1, nrow(Q_space)), 
      eps = 0
  )
  ## Change to same class as outputed by INLA::inla.scale.model
  Q_space <- as(Q_space, "dgTMatrix")   
}
