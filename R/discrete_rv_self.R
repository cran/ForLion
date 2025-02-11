#' function to generate discrete uniform random variables for initial random design points in ForLion
#'
#' @param n number of discrete random variables
#' @param xlist list of levels for variables to be generated
#'
#' @return list of discrete uniform random variables
#' @export
#'
#' @examples
#' n=3 #three discrete random variables
#' xlist=list(c(-1,1),c(-1,1),c(-1,0,1)) #two binary and one three-levels
#' discrete_rv_self(n, xlist)
#'
#'


discrete_rv_self = function(n, xlist){
  if(n==0){return(NA)}
  r.v = rep(NA, n)
  for(i in 1:n){
    r.v[i] = sample(unlist(xlist[i]), 1, replace = TRUE)
  }
  return(r.v)
}
