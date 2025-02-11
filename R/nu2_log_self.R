#' function to calculate the second derivative of nu function given eta for log link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return the second derivative of nu function for log link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu2_log_self(eta)
#'

nu2_log_self <- function(x) {
  exp(x);
}
