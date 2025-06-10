#' functions to solve 4th order polynomial function given coefficients
#'
#' @param c0 constant coefficient of polynomial function
#' @param c1 coefficient of 1st order term
#' @param c2 coefficient of 2nd order term
#' @param c3 coefficient of 3rd order term
#' @param c4 coefficient of 4th order term
#'
#' @return sol   the 4 solutions of the polynomial function
#' @export
#'
#' @examples polynomial_sol_J5(19,-53,19,-21,30)
polynomial_sol_J5 <- function(c0, c1, c2, c3, c4) {
  a0 = as.complex(c0/c4);
  a1 = as.complex(c1/c4);
  a2 = as.complex(c2/c4);
  a3 = as.complex(c3/c4);


  F1=27*a1^2-72*a0*a2+2*a2^3-9*a1*a2*a3+27*a0*a3^2

  E1=12*a0+a2^2-3*a1*a3

  D1=(F1+sqrt(F1^2-4*E1^3))^(1/3)

  G1=D1+(2^(2/3))*E1/D1

  A1=-2*a2/3+a3^2/4+G1/(3*2^(1/3))

  C1=-4*a2/3+(a3^2)/2-G1/(3*2^(1/3))+(-8*a1+4*a2*a3-a3^3)/(4*sqrt(A1))

  B1=-4*a2/3+(a3^2)/2-G1/(3*2^(1/3))-(-8*a1+4*a2*a3-a3^3)/(4*sqrt(A1))

  y11 = -(a3/4)-1/2*sqrt(A1) - 1/2*sqrt(B1);
  y12 = -(a3/4)-1/2*sqrt(A1) + 1/2*sqrt(B1);
  y13 = -(a3/4)+1/2*sqrt(A1) - 1/2*sqrt(C1);
  y14 = -(a3/4)+1/2*sqrt(A1) + 1/2*sqrt(C1);

  sol = c(y11, y12, y13, y14);
  sol;
}

