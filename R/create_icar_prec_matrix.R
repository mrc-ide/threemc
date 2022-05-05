#' @title Create Precision Matrix for ICAR Process
#'
#' @description Create the precision matrix for an ICAR process.
#'
#' @param sf_obj Shapefiles needed for adjacency, Default: NULL
#' @param area_lev  PSNU area level for specific country.
#' @param row.names Unique IDs for the areas, Default: NULL
#' @return ICAR precision matrix.
#'
#' @seealso
#'  \code{\link[spdep]{poly2nb}}
#'  \code{\link[spdep]{nb2mat}}
#'  \code{\link[naomi]{scale_gmrf_precision}}
#'  \code{\link[methods]{as}}
#' @rdname create_icar_prec_matrix
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
create_icar_prec_matrix <- function(sf_obj = NULL,
                                    area_lev = NULL,
                                    row.names = NULL) {
  if (is.null(area_lev)) {
    message(
      "area_lev missing, taken as maximum area level in sf_obj"
    )
    area_lev <- max(sf_obj$area_level, na.rm = TRUE)
  }
 
  sf_obj <- sf_obj %>%
    dplyr::filter(.data$area_level == area_lev)
  
  # if area_lev == 0, adjacency matrix will be a 1x1 matrix with single entry 0
  if (area_lev > 0) {
    # Creating neighbourhood structure
    Q_space <- spdep::poly2nb(sf_obj, row.names = sf_obj[, row.names])
    # Converting to adjacency matrix
    Q_space <- spdep::nb2mat(Q_space, style = "B", zero.policy = TRUE)
    
    # for precision matrix
    Q <- diag(rowSums(as.matrix(Q_space))) - 0.99 * Q_space
  } else {
    Q_space <- Matrix::Matrix(data = 0, nrow = 1, ncol = 1)
    Q <- as.matrix(0)
  }
  
  # Converting to sparse matrix
  Q_space <- methods::as(Q_space, "sparseMatrix")
  
  # Creating precision matrix from adjacency
  Q_space <- naomi::scale_gmrf_precision(
    Q   = Q, 
    A   = matrix(1, 1, nrow(Q_space)),
    eps = 0
  )
  
  # Change to same class as outputted by INLA::inla.scale.model
  Q_space <- methods::as(Q_space, "dgTMatrix")
}
