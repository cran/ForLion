#' functions to solve 3th order polynomial function given coefficients
#'
#' @param c0 constant coefficient of polynomial function
#' @param c1 coefficient of 1st order term
#' @param c2 coefficient of 2nd order term
#' @param c3 coefficient of 3rd order term
#'
#' @return sol   the 3 solutions of the polynomial function
#' @export
#'
#' @examples polynomial_sol_J4(0,9,6,1)
#'
polynomial_sol_J4<- function(c0, c1, c2, c3) {
  a0 = as.complex(c0/c3);
  a1 = as.complex(c1/c3);
  a2 = as.complex(c2/c3);

  q <- a1 / 3 - a2^2 / 9
  r <- (a1 * a2 - 3 * a0) / 6 - a2^3 / 27

  discriminant <- q^3 + r^2

  real_cuberoot <- function(x) {
    sign(x) * abs(x)^(1/3)
  }

  if (abs(discriminant) < .Machine$double.eps) {
    # q^3+r^2=0
    root_part <- real_cuberoot(Re(r))
    z1 <- 2 * root_part - a2 / 3
    z2 <- -root_part - a2 / 3
    return(c(z1, z2, z2))

  } else if (Re(discriminant) >= 0 && abs(Im(discriminant)) < .Machine$double.eps) {
    # q^3+r^2>0
    s1 <- real_cuberoot(Re(r + sqrt(discriminant)))
    s2 <- real_cuberoot(Re(r - sqrt(discriminant)))

  } else {
    # q^3+r^2<0
    s1 <- (r + sqrt(discriminant))^(1/3)
    s2 <- (r - sqrt(discriminant))^(1/3)
  }

  root1 <- s1 + s2 - a2 / 3
  root2 <- -(s1 + s2) / 2 - a2 / 3 + 1i * sqrt(3) * (s1 - s2) / 2
  root3 <- -(s1 + s2) / 2 - a2 / 3 - 1i * sqrt(3) * (s1 - s2) / 2

  sol <- c(root1, root2, root3)
  return(sol)
}

