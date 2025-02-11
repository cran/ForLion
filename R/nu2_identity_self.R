#' function to calculate the second derivative of nu function given eta for identity link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return the second derivative of nu function for identity link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu2_identity_self(eta)
#'

nu2_identity_self <- function(x) {
  rep(0, length(x))
}
