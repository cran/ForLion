#' function to calculate w = nu(eta) given eta for identity link
#'
#' @param x Numeric vector of eta, eta = X*beta.
#' @return A numeric vector representing the diagonal elements of the W matrix (nu(eta)).
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu_identity_self(eta)
#'

nu_identity_self <- function(x) {     # for identity link
  x;
}
