#' function to calculate w = nu(eta) given eta for logit link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return diagonal element of W matrix which is nu(eta)
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu_logit_self(eta)
#'
#'

nu_logit_self <- function(x) {  # ## functions to calculate w = nu(eta) given eta for logit link
  exp(x)/(1 + exp(x))^2;
}
