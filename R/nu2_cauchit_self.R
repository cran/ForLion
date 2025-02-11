#' function to calculate the second derivative of nu function given eta for cauchit link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return the second derivative of nu function for cauchit link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu2_cauchit_self(eta)
#'
#'

nu2_cauchit_self <- function(x) {
  16/(1+x^2)^4*((5*x^2-1)/(pi^2-4*atan(x)^2) - (20*x*atan(x)+6)/(pi^2-4*atan(x)^2)^2 + 8*pi^2/(pi^2-4*atan(x)^2)^3);
}
