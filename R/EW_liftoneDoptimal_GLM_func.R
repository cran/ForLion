#' EW Lift-one algorithm for D-optimal approximate design
#'
#' @param X Model matrix, with nrow = num of design points and ncol = num of parameters
#' @param E_w Diagonal of E_W matrix in Fisher information matrix, can be calculated EW_Xw_maineffects_self() function in the ForLion package
#' @param reltol reltol The relative convergence tolerance, default value 1e-5
#' @param maxit The maximum number of iterations, default value 100
#' @param random TRUE or FALSE, if TRUE then the function will run with additional "nram" number of initial allocation p00, default to be TRUE
#' @param nram When random == TRUE, the function will generate nram number of initial points, default is 3
#' @param p00 Specified initial design approximate allocation; default to be NULL, this will generate a random initial design
#'
#' @return p  EW D-optimal approximate allocation
#' @return p0 Initial approximate allocation that derived the reported EW D-optimal approximate allocation
#' @return Maximum The maximum of the determinant of the Fisher information matrix of the reported EW D-optimal design
#' @return convergence Convergence TRUE or FALSE
#' @return itmax number of the iteration
#' @export
#'
#' @examples
#' hfunc.temp = function(y) {c(y,1);};   # y -> h(y)=(y1,y2,y3,1)
#' link.temp="logit"
#' paras_lowerbound<-rep(-Inf, 4)
#' paras_upperbound<-rep(Inf, 4)
#' gjoint_b<- function(x) {
#' mu1 <- -0.5; sigma1 <- 1
#' mu2 <- 0.5; sigma2 <- 1
#' mu3 <- 1; sigma3 <- 1
#' mu0 <- 1; sigma0 <- 1
#' d1 <- stats::dnorm(x[1], mean = mu1, sd = sigma1)
#' d2 <- stats::dnorm(x[2], mean = mu2, sd = sigma2)
#' d3 <- stats::dnorm(x[3], mean = mu3, sd = sigma3)
#' d4 <- stats::dnorm(x[4], mean = mu0, sd = sigma0)
#' return(d1 * d2 * d3 * d4)
#' }
#' x.temp=matrix(data=c(-2,-1,-3,2,-1,-3,-2,1,-3,2,1,-3,-2,-1,3,2,-1,3,-2,1,3,2,1,3),ncol=3,byrow=TRUE)
#' m.temp=dim(x.temp)[1]     # number of design points
#' p.temp=length(paras_upperbound)    # number of predictors
#' Xmat.temp=matrix(0, m.temp, p.temp)
#' EW_wvec.temp=rep(0, m.temp)
#' for(i in 1:m.temp) {
#' htemp=EW_Xw_maineffects_self(x=x.temp[i,],joint_Func_b=gjoint_b, Lowerbounds=paras_lowerbound,
#'                              Upperbounds=paras_upperbound, link=link.temp, h.func=hfunc.temp);
#' Xmat.temp[i,]=htemp$X;
#' EW_wvec.temp[i]=htemp$E_w;
#' }
#' EW_liftoneDoptimal_GLM_func(X=Xmat.temp, E_w=EW_wvec.temp, reltol=1e-8, maxit=1000,
#'                             random=TRUE, nram=3, p00=NULL)

EW_liftoneDoptimal_GLM_func <- function(X, E_w, reltol=1e-5, maxit=100, random=FALSE, nram=3, p00=NULL)  {
  ## E_w=E_w[1,2,...,m] are strictly positive
  # if random=TRUE, run 5 random initial points and pick up the best; default initial p1=p2=...=1/m
  # output: p=p--optimal design based on "det"
  #         Maximum--maximized value of "det"
  #         convergence -- "TRUE" indicates success
  #         p0 -- initial p
  #         itmax -- number of iterations
  m = dim(X)[1];
  d = dim(X)[2];
  E_w = E_w[1:m];
  if(min(E_w) <= 0) {
    warning("\nE_w's need to be strictly positive!\n");
    return(0);
  };
  ftemp <- function(p) { det(t(X * (p*E_w)) %*% X);};
  if(is.null(p00)) p00=rep(1/m,m);
  maximum = ftemp(p00);
  maxvec = stats::rexp(m);
  convergence = FALSE;
  p = p00;
  ind = 0;
  while((ind < maxit) && ((max(maxvec,na.rm=TRUE)/min(maxvec,na.rm=TRUE))-1 > reltol)) {
    io = sample(x=(1:m), size=m);
    for(i in 1:m) {    # run updating in random order of E_w
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
  p.ans=p; maximum.ans=maximum; if((max(maxvec,na.rm=TRUE)/min(maxvec,na.rm=TRUE))-1 <= reltol) convergence=TRUE;itmax=ind;
  if(random) for(j in 1:nram) {
    p0=stats::rexp(m);p0=p0/sum(p0);
    p=p0;
    maxvec = stats::rexp(m);
    maximum = ftemp(p);
    ind = 0;
    while((ind < maxit) && ((max(maxvec,na.rm=TRUE)/min(maxvec,na.rm=TRUE))-1 > reltol)) {
      io = sample(x=(1:m), size=m);
      for(i in 1:m) {    # run updating in random order of E_w
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
      convergence=FALSE;
      if((max(maxvec,na.rm=TRUE)/min(maxvec,na.rm=TRUE))-1 <= reltol) convergence=TRUE;
      p00=p0;
      itmax=ind;
    }
  }
  # list(p=p.ans, p0=p00, Maximum=maximum.ans, convergence=convergence, itmax=itmax);  # convergence=T indicates success
  output<-list(p=p.ans, p0=p00, Maximum=maximum.ans, convergence=convergence, itmax=itmax);
  class(output) <- "list_output"
  return(output)
  }
