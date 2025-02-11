#' function to calculate w = nu(eta) given eta for loglog link
#'
#' @param x vector of eta, eta=X*beta
#'
#' @return diagonal element of W matrix which is nu(eta)
#' @export
#'
#' @examples
#' eta = c(1,2,3,4)
#' nu_loglog_self(eta)
#'
#'

nu_loglog_self <- function(x) {  # for c-log-log and log-log links
  exp(2*x)/(exp(exp(x))-1);
}
