#' function to calculate first derivative of nu function given eta for identity link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return the first derivative of nu function given eta for identity link
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu1_identity_self(eta)
#'

nu1_identity_self <- function(x) {
  rep(1, length(x)) #constant
}
