#' functions to solve 2th order polynomial function given coefficients
#'
#' @param c0 constant coefficient of polynomial function
#' @param c1 coefficient of 1st order term
#' @param c2 coefficient of 2nd order term
#'
#' @return sol   the 2 solutions of the polynomial function
#' @export
#'
#' @examples polynomial_sol_J3(-2,-3, 1)
#'
polynomial_sol_J3<- function(c0, c1, c2) {
  c2 <- as.complex(c2)
  c1 <- as.complex(c1)
  c0 <- as.complex(c0)

  discriminant <- c1^2 - 4*c2*c0

  root1 <- (-c1 + sqrt(discriminant))/(2*c2)
  root2 <- (-c1 - sqrt(discriminant))/(2*c2)

  sol <- c(root1, root2)
  return(sol)
}


