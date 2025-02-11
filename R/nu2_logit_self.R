#' function to calculate the second derivative of nu function given eta for logit link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return the second derivative of nu function for logit link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu2_logit_self(eta)
#'

nu2_logit_self <- function(x) {
  exp(x)*(1-4*exp(x)+exp(2*x))/(1 + exp(x))^4;
}
