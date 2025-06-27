#' Lift-one algorithm for D-optimal approximate design
#' @importFrom stats rexp
#' @param X Model matrix, with nrow = num of design points and ncol = num of parameters
#' @param w Diagonal of W matrix in Fisher information matrix, can be calculated Xw_maineffects_self() function in the ForLion package
#' @param reltol The relative convergence tolerance, default value 1e-5
#' @param maxit The maximum number of iterations, default value 100
#' @param random TRUE or FALSE, if TRUE then the function will run lift-one with additional "nram" number of random approximate allocation, default to be FALSE
#' @param nram when random == TRUE, the function will run lift-one nram number of initial proportion p00, default is 3
#' @param p00 Specified initial design approximate allocation; default to be NULL, this will generate a random initial design
#'
#' @return p              D-optimal approximate allocation
#' @return p0             Initial approximate allocation that derived the reported D-optimal approximate allocation
#' @return Maximum        The maximum of the determinant of the Fisher information matrix of the reported D-optimal design
#' @return convergence    Convergence TRUE or FALSE
#' @return itmax          number of the iteration
#' @export
#'
#' @examples
#' hfunc.temp = function(y) {c(y,y[4]*y[5],1);};   # y -> h(y)=(y1,y2,y3,y4,y5,y4*y5,1)
#' link.temp="logit"
#' x.temp = matrix(data=c(25.00000,1,-1,1,-1,25.00000,1,1,1,-1,32.06741,-1,1,-1,1,40.85698,
#' -1,1,1,-1,28.86602,-1,1,-1,-1,29.21486,-1,-1,1,1,25.00000,1,1,1,1, 25.00000,1,1,-1,-1),
#' ncol=5, byrow=TRUE)
#' b.temp = c(0.3197169,  1.9740922, -0.1191797, -0.2518067,  0.1970956,  0.3981632, -7.6648090)
#' X.mat = matrix(,nrow=8, ncol=7)
#' w.vec = rep(NA,8)
#' for(i in 1:8) {
#' htemp=Xw_maineffects_self(x=x.temp[i,], b=b.temp, link=link.temp, h.func=hfunc.temp);
#' X.mat[i,]=htemp$X;
#' w.vec[i]=htemp$w;
#' };
#' liftoneDoptimal_GLM_func(X=X.mat, w=w.vec, reltol=1e-5, maxit=500, random=TRUE, nram=3, p00=NULL)
#'



liftoneDoptimal_GLM_func  <- function(X, w, reltol=1e-5, maxit=100, random=FALSE, nram=3, p00=NULL)  {
  ## Lift-one algorithm for D-optimal approximate design step(iii) and after
  ## Version: 09/29/2018
  ## Find approximate optimal solution p maximizing |X'WX| given X and w
  ## X[m*d] = X[A1,...,Ad]
  ## p00 can be specified
  ## w=w[1,2,...,m] are strictly positive
  # if random=T, run 5 random initial points and pick up the best; default initial p1=p2=...=1/m
  # output: p=p--optimal design based on "det"
  #         Maximum--maximized value of "det"
  #         convergence -- "T" indicates success
  #         p0 -- initial p
  #         itmax -- number of iterations
  m = dim(X)[1];
  d = dim(X)[2];
  w = w[1:m];
  if(min(w) <= 0) {
    warning("\nW's need to be strictly positive!\n");
    return(0);
  };
  ftemp <- function(p) { det(t(X * (p*w)) %*% X);};
  if(is.null(p00)) p00=rep(1/m,m);
  maximum = ftemp(p00);
  maxvec = stats::rexp(m);
  convergence = F;
  p = p00;
  ind = 0;
  while((ind < maxit) && ((max(maxvec,na.rm=T)/min(maxvec,na.rm=T))-1 > reltol)) {
    io = sample(x=(1:m), size=m);
    for(i in 1:m) {    # run updating in random order of w
      if(p[io[i]]>0) {
        ptemp1 = p/(1-p[io[i]]);
        ptemp1[io[i]] = 0;
        b = ftemp(ptemp1);      # b=fs(0)
        a = (maximum - b*(1-p[io[i]])^d)/(p[io[i]]*(1-p[io[i]])^(d-1));
      } else {         # p[io[i]]==0
        b = maximum;
        ptemp1 = p/2;
        ptemp1[io[i]] = 1/2;          # for fs(1/2)
        a = ftemp(ptemp1)*2^d - b;
      }
      if(a > b*d) x=(a-b*d)/((a-b)*d) else x=0;
      ptemp1 = p*(1-x)/(1-p[io[i]]);
      ptemp1[io[i]] = x;
      if(a > b*d) maximum = ((d-1)/(a-b))^(d-1)*(a/d)^d else maximum=b;
      p = ptemp1;
      maxvec[io[i]] = maximum;
    }
    ind = ind+1;
  }
  p.ans=p; maximum.ans=maximum; if((max(maxvec,na.rm=T)/min(maxvec,na.rm=T))-1 <= reltol) convergence=T;itmax=ind;
  if(random) for(j in 1:nram) {
    p0=stats::rexp(m);p0=p0/sum(p0);
    p=p0;
    maxvec = stats::rexp(m);
    maximum = ftemp(p);
    ind = 0;
    while((ind < maxit) && ((max(maxvec,na.rm=T)/min(maxvec,na.rm=T))-1 > reltol)) {
      io = sample(x=(1:m), size=m);
      for(i in 1:m) {    # run updating in random order of w
        if(p[io[i]]>0) {
          ptemp1 = p/(1-p[io[i]]);
          ptemp1[io[i]] = 0;
          b = ftemp(ptemp1);      # b=fs(0)
          a = (maximum - b*(1-p[io[i]])^d)/(p[io[i]]*(1-p[io[i]])^(d-1));
        } else {         # p[io[i]]==0
          b = maximum;
          ptemp1 = p/2;
          ptemp1[io[i]] = 1/2;          # for fs(1/2)
          a = ftemp(ptemp1)*2^d - b;
        }
        if(a > b*d) x=(a-b*d)/((a-b)*d) else x=0;
        ptemp1 = p*(1-x)/(1-p[io[i]]);
        ptemp1[io[i]] = x;
        if(a > b*d) maximum = ((d-1)/(a-b))^(d-1)*(a/d)^d else maximum=b;
        p = ptemp1;
        maxvec[io[i]] = maximum;
      }
      ind = ind+1;
    }
    if(maximum > maximum.ans) {
      maximum.ans=maximum;
      p.ans=p;
      convergence=F;
      if((max(maxvec,na.rm=T)/min(maxvec,na.rm=T))-1 <= reltol) convergence=T;
      p00=p0;
      itmax=ind;
    }
  }
  # list(p=p.ans, p0=p00, Maximum=maximum.ans, convergence=convergence, itmax=itmax);  # convergence=T indicates success

  #define S3 class
  output<-list(p=p.ans, p0=p00, Maximum=maximum.ans, convergence=convergence, itmax=itmax);
  class(output) <- "list_output"
  return(output)
}
