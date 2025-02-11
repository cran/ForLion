#' function to calculate the first derivative of nu function given eta for probit link
#' @param x vector of eta, eta=X*beta
#'
#' @return the first derivative of nu function for probit link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu1_probit_self(eta)
#'


nu1_probit_self <- function(x) {
  -2*x*stats::dnorm(x)^2/(stats::pnorm(x)*(1-stats::pnorm(x))) + stats::dnorm(x)^3*(2*stats::pnorm(x)-1)/((1-stats::pnorm(x))^2*stats::pnorm(x)^2);
}
