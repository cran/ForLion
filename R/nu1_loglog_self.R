#' function to calculate the first derivative of nu function given eta for log-log link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return the first derivative of nu function given eta for log-log link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu1_loglog_self(eta)
#'


nu1_loglog_self <- function(x) {
  2*exp(2*x)/(exp(exp(x))-1) - exp(3*x+exp(x))/(exp(exp(x))-1)^2;
}
