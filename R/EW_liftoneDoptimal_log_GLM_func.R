#' EW Lift-one algorithm for D-optimal approximate design in log scale
#'
#' @param X Model matrix, with nrow = num of design points and ncol = num of parameters
#' @param E_w Diagonal of E_W matrix in Fisher information matrix, can be calculated EW_Xw_maineffects_self() function in the ForLion package
#' @param reltol reltol The relative convergence tolerance, default value 1e-5
#' @param maxit The maximum number of iterations, default value 100
#' @param random TRUE or FALSE, if TRUE then the function will run lift-one with additional "nram" number of random approximate allocation, default to be FALSE
#' @param nram when random == TRUE, the function will run lift-one nram number of initial proportion p00, default is 3
#' @param p00 Specified initial design approximate allocation; default to be NULL, this will generate a random initial design
#'
#' @return p             EW D-optimal approximate allocation
#' @return p0            Initial approximate allocation that derived the reported EW D-optimal approximate allocation
#' @return Maximum       The maximum of the determinant of the Fisher information matrix of the reported EW D-optimal design
#' @return convergence   Convergence TRUE or FALSE
#' @return itmax         number of the iteration
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
#' x.temp=matrix(data=c(-2,-1,-3,2,-1,-3,-2,1,-3,2,1,-3,-2,-1,3,2,-1,3,-2,1,3,2,1,3),
#'               ncol=3,byrow=TRUE)
#' m.temp=dim(x.temp)[1]     # number of design points
#' p.temp=length(paras_upperbound)    # number of predictors
#' Xmat.temp=matrix(0, m.temp, p.temp)
#' EW_wvec.temp=rep(0, m.temp)
#' for(i in 1:m.temp) {
#' htemp=EW_Xw_maineffects_self(x=x.temp[i,],Integral_based=TRUE,joint_Func_b=gjoint_b,
#' Lowerbounds=paras_lowerbound, Upperbounds=paras_upperbound, link=link.temp,
#' h.func=hfunc.temp);
#' Xmat.temp[i,]=htemp$X;
#' EW_wvec.temp[i]=htemp$E_w;
#' }
#' EW_liftoneDoptimal_GLM_func(X=Xmat.temp, E_w=EW_wvec.temp, reltol=1e-8, maxit=1000, random=TRUE,
#'                             nram=3, p00=NULL)

EW_liftoneDoptimal_log_GLM_func <- function(X, E_w, reltol=1e-5, maxit=100, random=FALSE, nram=3, p00=NULL)  {
  ## E_w=E_w[1,2,...,m] are strictly positive
  ldiff <- function(logx, logy) {
    ifelse(logx > logy, log(1 - exp(logy-logx)) + logx, -Inf);
  }
  m = dim(X)[1];
  d = dim(X)[2];
  E_w = E_w[1:m];
  if(min(E_w) <= 0) {
    warning("\nE_w's need to be strictly positive!\n");
    return(0);
  };
  ftemp <- function(p) { determinant(t(X * (p*E_w)) %*% X)$modulus[1];};  # log det
  if(is.null(p00)) p00=rep(1/m,m);
  maximum = ftemp(p00);
  maxvec = stats::rexp(m);
  convergence = FALSE;
  p = p00;
  ind = 0;
  while((ind < maxit) && (max(maxvec)-min(maxvec) > reltol)) {
    io = sample(x=(1:m), size=m);
    for(i in 1:m) {    # run updating in random order of E_w
      if(p[io[i]]>0) {
        ptemp1 = p;
        ptemp1[io[i]] = 0;
        b1 = ftemp(ptemp1);      # b1=fs(0)*(1-p_i)^d=b*fs(0)*(1-p_i)^d, log scale
        if(maximum <= b1) x=0 else {
          a = ldiff(maximum, b1) - log(p[io[i]]) - (d-1)*log(1-p[io[i]]); # in log scale
          b = ftemp(ptemp1/(1-p[io[i]]));
          if(a <= b + log(d)) x=0 else {
            x=exp(ldiff(a, b+log(d)) - ldiff(a, b) - log(d));
          };
        };
      } else {         # p[io[i]]==0
        b = maximum;   # log scale
        ptemp1 = p;
        ptemp1[io[i]] = 1;
        a1 = ftemp(ptemp1);     # a1=fs(1/2)*2^d, log scale
        if(a1 <= b) x=0 else {
          a = ldiff(a1, b);
          if(a <= b + log(d)) x=0 else {
            x=exp(ldiff(a, b+log(d)) - ldiff(a, b) - log(d));
          };
        };
      };               # end of p[io[i]]=0
      ptemp1 = p*(1-x)/(1-p[io[i]]);
      ptemp1[io[i]] = x;
      p = ptemp1;
      maximum=ftemp(p);
      maxvec[io[i]] = maximum;
    };                 # end of "for(i in 1:m)"
    ind = ind+1;
    #delete after
    #cat("\n", ind,"iteration: \n p=",p, "\n maximum=", maximum, "\n x=", x)
  }

  p.ans=p; maximum.ans=maximum;
  if(max(maxvec)-min(maxvec) <= reltol) convergence=TRUE;
  itmax=ind;
  if(random) for(j in 1:nram) {
    p0=stats::rexp(m);p0=p0/sum(p0);
    p=p0;
    maxvec = stats::rexp(m);
    maximum = ftemp(p);
    ind = 0;
    while((ind < maxit) && (max(maxvec)-min(maxvec) > reltol)) {
      io = sample(x=(1:m), size=m);
      for(i in 1:m) {    # run updating in random order of E_w
        if(p[io[i]]>0) {
          ptemp1 = p;
          ptemp1[io[i]] = 0;
          b1 = ftemp(ptemp1);      # b1=fs(0)*(1-p_i)^d=b*fs(0)*(1-p_i)^d, log scale
          if(maximum <= b1) x=0 else {
            a = ldiff(maximum, b1) - log(p[io[i]]) - (d-1)*log(1-p[io[i]]); # in log scale
            b = ftemp(ptemp1/(1-p[io[i]]));
            if(a <= b + log(d)) x=0 else {
              x=exp(ldiff(a, b+log(d)) - ldiff(a, b) - log(d));
            };
          };
        } else {         # p[io[i]]==0
          b = maximum;   # log scale
          ptemp1 = p;
          ptemp1[io[i]] = 1;
          a1 = ftemp(ptemp1);     # a1=fs(1/2)*2^d, log scale
          if(a1 <= b) x=0 else {
            a = ldiff(a1, b);
            if(a <= b + log(d)) x=0 else {
              x=exp(ldiff(a, b+log(d)) - ldiff(a, b) - log(d));
            };
          };
        };               # end of p[io[i]]=0
        ptemp1 = p*(1-x)/(1-p[io[i]]);
        ptemp1[io[i]] = x;
        p = ptemp1;
        maximum=ftemp(p);
        maxvec[io[i]] = maximum;
      };                 # end of "for(i in 1:m)"
      ind = ind+1;
    }
    if(maximum > maximum.ans) {
      maximum.ans=maximum;
      p.ans=p;
      convergence=FALSE;
      if(max(maxvec)-min(maxvec) <= reltol) convergence=TRUE;
      p00=p0;
      itmax=ind;
    }
  }
  #list(p=p.ans, p0=p00, Maximum=exp(maximum.ans), convergence=convergence, itmax=itmax);  # convergence=T indicates success

  #define S3 class
  output<-list(p=p.ans, p0=p00, Maximum=exp(maximum.ans), convergence=convergence, itmax=itmax);
  class(output) <- "list_output"
  return(output)
  }
