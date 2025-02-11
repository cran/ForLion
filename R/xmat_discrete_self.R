#' Generate GLM random initial designs within ForLion algorithm
#'
#' @param xlist a list of factor levels within ForLion algorithm, for example, a binary factor might be c(-1,1), a continuous factor within range of (25,45) will be c(25, 45).
#' @param rowmax maximum number of rows of the design matrix
#'
#' @return design matrix of all possible combinations of discrete factors levels with min and max of the continuous factors.
#' @export
#'
#' @examples
#' #define list of factor levels for one continuous factor, four binary factors
#' factor.level.temp = list(c(25,45),c(-1,1),c(-1,1),c(-1,1),c(-1,1))
#' xmat_discrete_self(xlist = factor.level.temp)
#'
#'


xmat_discrete_self <- function(xlist, rowmax=NULL) {
  dx=length(xlist);   # number of factors
  mvec=rep(0,dx);     # (m1,m2,...,md)
  for(i in 1:dx) mvec[i]=length(xlist[[i]]);
  nx=prod(mvec);      # total number of level combinations
  xmat=matrix(0, nx, dx);   # matrix of level combinations for output
  xtemp=(0:(nx-1));
  ntemp=nx;
  for(i in 1:dx) {
    ntemp=ntemp/mvec[i];
    xmat[,i]=floor(xtemp/ntemp);
    xtemp=xtemp-xmat[,i]*ntemp;
    xmat[,i]=xlist[[i]][xmat[,i]+1];
  };
  if(!is.null(rowmax)) {
    rtemp=min(rowmax, nx);
    itemp=sort(sample(1:nx, size=rtemp, replace=FALSE));
    xmat=xmat[itemp,];
  };
  xmat;
}
