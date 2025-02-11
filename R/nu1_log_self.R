#' function to calculate first derivative of nu function given eta for log link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return the first derivative of nu function given eta for log link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu1_log_self(eta)
#'

nu1_log_self <- function(x) {
  exp(x);
}
