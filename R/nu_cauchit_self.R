#' function to calculate w = nu(eta) given eta for cauchit link
#'
#' @param x a list of eta = X*beta
#'
#' @return diagonal element of W matrix which is nu(eta)
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu_cauchit_self(eta)
#'
#'
nu_cauchit_self <- function(x) { # for cauchit link
  1/(1+x^2)^2/(pi^2/4-atan(x)^2);
}
