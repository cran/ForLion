#' function to calculate the second derivative of nu function given eta for probit link
#' @importFrom stats dnorm pnorm
#' @param x vector of eta, eta=X*beta
#'
#' @return the second derivative of nu function for probit link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu2_probit_self(eta)
#'


nu2_probit_self <- function(x) {
  2*(2*x^2-1)*stats::dnorm(x)^2/(stats::pnorm(x)*(1-stats::pnorm(x))) - (5*x*stats::dnorm(x)^3*(2*stats::pnorm(x)-1)+6*stats::dnorm(x)^4)/((1-stats::pnorm(x))^2*stats::pnorm(x)^2) + 2*stats::dnorm(x)^4/((1-stats::pnorm(x))^3*stats::pnorm(x)^3);
}
