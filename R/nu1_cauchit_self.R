#' function to calculate first derivative of nu function given eta for cauchit link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return the first derivative of nu function given eta for cauchit link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu1_cauchit_self(eta)
#'

nu1_cauchit_self <- function(x) {
  16/(1+x^2)^3*(-x/(pi^2-4*atan(x)^2) + 2*atan(x)/(pi^2-4*atan(x)^2)^2);
}
