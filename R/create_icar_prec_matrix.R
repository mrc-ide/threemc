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
#'  \code{\link[INLA]{inla.scale.model}}
#' @rdname create_icar_prec_matrix
#' @export
#' 
#' @importFrom spdep poly2nb nb2mat
#' @importFrom INLA inla.scale.model
create_icar_prec_matrix <- function(sf_obj = NULL,
                                    row.names = NULL) {
  # Creating neighbourhood structure
  Q_space <-  poly2nb(sf_obj, row.names = sf_obj[, row.names])
  # Converting to adjacency matrix
  Q_space <- nb2mat(Q_space, style = "B", zero.policy = TRUE)
  # Converting to sparse matrix
  Q_space <- as(Q_space, "sparseMatrix")
  # Creating precision matrix from adjacency
  Q_space <- inla.scale.model(
    diag(rowSums(Q_space)) - 0.99 * Q_space,
    constr = list(A = matrix(1, 1, nrow(Q_space)), e = 0)
  )
}
