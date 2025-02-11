#' function to calculate the first derivative of nu function given eta for logit link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return the first derivative of nu function given eta for logit link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu1_logit_self(eta)
#'
#'


nu1_logit_self <- function(x) {
  -exp(x)*(exp(x)-1)/(1 + exp(x))^3;
}
