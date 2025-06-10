#' SVD Inverse Of A Square Matrix
#' @description
#' This function returns the inverse of a matrix using singular value decomposition.
#' If the matrix is a square matrix, this should be equivalent to using the solve function.
#' If the matrix is not a square matrix, then the result is the Moore-Penrose pseudo inverse.
#'
#' @param x the matrix for calculation of inverse
#'
#' @return the inverse of the matrix x
#' @export
#'
#' @examples
#' x = diag(4)
#' svd_inverse(x)
#'



svd_inverse <- function (x) #from R package "matrixcalc"
{
  if (!is.numeric(x)) {
    stop("argument x is not numeric")
  }
  if (is.matrix(x)) {
    Xmat <- x
  }
  else {
    if (is.vector(x)) {
      Xmat <- matrix(x, nrow = length(x), ncol = 1)
    }
    else {
      stop("argument x is neither a vector nor a matrix")
    }
  }
  svdXmat <- svd(Xmat)
  U <- svdXmat$u
  d <- svdXmat$d
  if (any(d == 0)) {
    stop("x has at least one zero singular value")
  }
  if (length(d) == 1) {
    Dinv <- matrix(1/d, nrow = 1)
  }
  else {
    Dinv <- diag(1/d)
  }
  V <- svdXmat$v
  Xinv <- V %*% Dinv %*% t(U)
  return(Xinv)
}
