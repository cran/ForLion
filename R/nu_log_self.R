#' function to calculate w = nu(eta) given eta for log link
#'
#' @param x Numeric vector of eta, eta = X*beta.
#' @return A numeric vector representing the diagonal elements of the W matrix (nu(eta)).
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu_log_self(eta)
#'

nu_log_self <- function(x) {     # for log link
  exp(x);
}
