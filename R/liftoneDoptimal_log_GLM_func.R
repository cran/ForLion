#' Lift-one algorithm for D-optimal approximate design in log scale
#' @importFrom stats rexp
#' @param X Model matrix, with nrow = num of design points and ncol = num of parameters
#' @param w Diagonal of W matrix in Fisher information matrix, can be calculated Xw_maineffects_self() function in the ForLion package
#' @param reltol The relative convergence tolerance, default value 1e-5
#' @param maxit The maximum number of iterations, default value 100
#' @param random TRUE or FALSE, if TRUE then the function will run with additional "nram" number of initial allocation p00, default to be TRUE
#' @param nram When random == TRUE, the function will generate nram number of initial points, default is 3
#' @param p00 Specified initial design approximate allocation; default to be NULL, this will generate a random initial design
#'
#' @return p               D-optimal approximate allocation
#' @return p0              Initial approximate allocation that derived the reported D-optimal approximate allocation
#' @return Maximum The     maximum of the determinant of the expected Fisher information matrix of the reported D-optimla design
#' @return convergence     Convergence TRUE or FALSE
#' @return itmax           number of the iteration
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
#' liftoneDoptimal_log_GLM_func(X=X.mat, w=w.vec, reltol=1e-5, maxit=500,
#' random=TRUE, nram=3, p00=NULL)




liftoneDoptimal_log_GLM_func  <- function(X, w, reltol=1e-5, maxit=100, random=FALSE, nram=3, p00=NULL)  {
  ldiff <- function(logx, logy) {
    ifelse(logx > logy, log(1 - exp(logy-logx)) + logx, -Inf);
  }
  m = dim(X)[1];
  d = dim(X)[2];
  w = w[1:m];
  if(min(w) <= 0) {
    warning("\nW's need to be strictly positive!\n");
    return(0);
  };
  ftemp <- function(p) { determinant(t(X * (p*w)) %*% X)$modulus[1];};  # log det
  if(is.null(p00)) p00=rep(1/m,m);
  maximum = ftemp(p00);
  maxvec = stats::rexp(m);
  convergence = F;
  p = p00;
  ind = 0;
  while((ind < maxit) && (max(maxvec)-min(maxvec) > reltol)) {
    io = sample(x=(1:m), size=m);
    for(i in 1:m) {    # run updating in random order of w
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
  if(max(maxvec)-min(maxvec) <= reltol) convergence=T;
  itmax=ind;
  if(random) for(j in 1:nram) {
    p0=stats::rexp(m);p0=p0/sum(p0);
    p=p0;
    maxvec = stats::rexp(m);
    maximum = ftemp(p);
    ind = 0;
    while((ind < maxit) && (max(maxvec)-min(maxvec) > reltol)) {
      io = sample(x=(1:m), size=m);
      for(i in 1:m) {    # run updating in random order of w
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
      convergence=F;
      if(max(maxvec)-min(maxvec) <= reltol) convergence=T;
      p00=p0;
      itmax=ind;
    }
  }
  # list(p=p.ans, p0=p00, Maximum=exp(maximum.ans), convergence=convergence, itmax=itmax);  # convergence=T indicates success

  #define S3 class
  output<-list(p=p.ans, p0=p00, Maximum=exp(maximum.ans), convergence=convergence, itmax=itmax);
  class(output) <- "list_output"
  return(output)
}

