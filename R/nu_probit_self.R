#' function to calculate w = nu(eta) given eta for probit link
#' @importFrom stats dnorm pnorm
#' @param x vector of eta, eta=X*beta
#'
#' @return diagonal element of W matrix which is nu(eta)
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu_probit_self(eta)
#'
#'

nu_probit_self <- function(x) {  # for probit link
  stats::dnorm(x)^2/(stats::pnorm(x)*(1-stats::pnorm(x)));
}
