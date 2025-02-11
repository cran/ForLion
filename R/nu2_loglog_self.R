#' function to calculate the second derivative of nu function given eta for loglog link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return the second derivative of nu function for loglog link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu2_loglog_self(eta)
#'


nu2_loglog_self <- function(x) {  # for c-log-log and log-log links
  exp(2*x)*(exp(x)-1)*(exp(x)-4)/(exp(exp(x))-1) + exp(3*x)*(3*exp(x)-5)/(exp(exp(x))-1)^2 + 2*exp(4*x)/(exp(exp(x))-1)^3;
}
