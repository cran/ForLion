#' ForLion for generalized linear models
#' @description
#' ForLion algorithm to find D-optimal design for GLM models with mixed factors, reference: Section 4 in Huang, Li, Mandal, Yang (2024).
#' Factors may include discrete factors with finite number of distinct levels and continuous factors with specified interval range (min, max), continuous factors, if any, must serve as main-effects only, allowing merging points that are close enough.
#' Continuous factors first then discrete factors, model parameters should in the same order of factors.
#' @param n.factor vector of numbers of distinct levels, "0" indicates continuous factors, "0"s always come first, "2" or above indicates discrete factor, "1" is not allowed
#' @param factor.level list of distinct levels, (min, max) for continuous factor, continuous factors first, should be the same order as n.factor
#' @param hfunc function for obtaining model matrix h(y) for given design point y, y has to follow the same order as n.factor
#' @param bvec assumed parameter values of model parameters beta, same length of h(y)
#' @param link link function, default "logit", other links: "probit", "cloglog", "loglog", "cauchit", "log", "identity"
#' @param reltol the relative convergence tolerance, default value 1e-5
#' @param rel.diff points with distance less than that will be merged, default value 0
#' @param maxit the maximum number of iterations, default value 100
#' @param random TRUE or FALSE, if TRUE then the function will run lift-one with additional "nram" number of random approximate allocation, default to be FALSE
#' @param nram when random == TRUE, the function will run lift-one nram number of initial proportion p00, default is 3
#' @param logscale TRUE or FALSE, if TRUE then the ForLion will run lift-one with logscale, which is liftoneDoptimal_log_GLM_func(); if FALSE then ForLion will run lift-one without logscale, which is liftoneDoptimal_GLM_func()
#' @param rowmax maximum number of points in the initial design, default NULL indicates no restriction
#' @param Xini initial list of design points, default NULL will generate random initial design points
#'
#' @return m number of design points
#' @return x.factor matrix with rows indicating design point
#' @return p D-optimal approximate allocation
#' @return det Optimal determinant of Fisher information matrix
#' @return convergence TRUE or FALSE
#' @return min.diff the minimum Euclidean distance between design points
#' @return x.close  a pair of design points with minimum distance
#' @return itmax iteration of the algorithm
#' @export
#'
#' @examples
#' #Example 3 in Huang, Li, Mandal, Yang (2024), electrostatic discharge experiment
#' hfunc.temp = function(y) {c(y,y[4]*y[5],1);};   # y -> h(y)=(y1,y2,y3,y4,y5,y4*y5,1)
#' n.factor.temp = c(0, 2, 2, 2, 2)  # 1 continuous factor with 4 discrete factors
#' factor.level.temp = list(c(25,45),c(-1,1),c(-1,1),c(-1,1),c(-1,1))
#' link.temp="logit"
#' b.temp = c(0.3197169,  1.9740922, -0.1191797, -0.2518067,  0.1970956,  0.3981632, -7.6648090)
#' ForLion_GLM_Optimal(n.factor=n.factor.temp, factor.level=factor.level.temp, hfunc=hfunc.temp,
#' bvec=b.temp, link=link.temp, reltol=1e-2, rel.diff=0.03, maxit=500, random=FALSE,
#' nram=3, logscale=TRUE)
#'


ForLion_GLM_Optimal <- function(n.factor, factor.level, hfunc, bvec, link, reltol=1e-5, rel.diff=0, maxit=100, random=FALSE, nram=3, logscale=FALSE, rowmax=NULL, Xini=NULL) {
  d.factor=length(n.factor);             # number of factors
  p.factor=length(bvec);                 # number of predictors
  k.continuous=sum(n.factor==0);         # number of continuous factors
  if(rel.diff==0) rel.diff=reltol;
  # functions for nu(eta), nu'(eta), nu''(eta) given eta
  nutemp=nu_logit_self; nu1temp=nu1_logit_self; nu2temp=nu2_logit_self;
  if(link=="probit") {nutemp=nu_probit_self; nu1temp=nu1_probit_self; nu2temp=nu2_probit_self;};
  if(link=="cloglog") {nutemp=nu_loglog_self; nu1temp=nu1_loglog_self; nu2temp=nu2_loglog_self;};
  if(link=="loglog") {nutemp=nu_loglog_self; nu1temp=nu1_loglog_self; nu2temp=nu2_loglog_self;};
  if(link=="cauchit") {nutemp=nu_cauchit_self; nu1temp=nu1_cauchit_self; nu2temp=nu2_cauchit_self;};
  if(link=="log") {nutemp=nu_log_self; nu1temp=nu1_log_self; nu2temp=nu2_log_self;};
  if(link=="identity") {nutemp=nu_identity_self; nu1temp=nu1_identity_self; nu2temp=nu2_identity_self;};
  #    Case I: all factors are discrete
  if(k.continuous==0) {
    if(is.null(Xini)) xtemp=xmat_discrete_self(factor.level, rowmax=rowmax) else xtemp=Xini;
    m.design=nrow(xtemp);                   # initial number of design points
    X.mat = matrix(0, m.design, p.factor);  # initial model matrix X
    w.vec = rep(0, m.design);     # w vector
    for(i in 1:m.design) {
      htemp=Xw_maineffects_self(x=xtemp[i,], b=bvec, link=link, h.func=hfunc);
      X.mat[i,]=htemp$X;
      w.vec[i]=htemp$w;
    };
    if(logscale) optemp=liftoneDoptimal_log_GLM_func (X=X.mat, w=w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram) else {
      optemp=liftoneDoptimal_GLM_func (X=X.mat, w=w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram);
    };
    m.design=sum(optemp$p>0);       # updated number of design point
    x.design=xtemp[optemp$p>0,];  # updated list of design points
    p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
    det.design=optemp$Maximum;      # optimized |X'WX|
    X.mat = X.mat[optemp$p>0,];   # updated model matrix
    w.vec = w.vec[optemp$p>0];    # updated w vector
    converge.design=optemp$convergence;  # TRUE or FALSE
    itmax.design=0;               # no candidate design points considered
  };                        #     End of Case I
  #    Case II: all factors are continuous
  if(k.continuous==d.factor) {
    lvec=uvec=rep(0, d.factor);     # lower bounds and upper bounds for continuous factors
    for(i in 1:d.factor) {lvec[i]=min(factor.level[[i]]); uvec[i]=max(factor.level[[i]]);};
    if(is.null(Xini)) xtemp=xmat_discrete_self(factor.level, rowmax=rowmax) else xtemp=Xini;
    if(k.continuous==1) m.design=length(xtemp) else m.design=nrow(xtemp);                   # initial number of design points
    X.mat = matrix(0, m.design, p.factor);  # initial model matrix X
    w.vec = rep(0, m.design);     # w vector
    for(i in 1:m.design) {
      if(k.continuous==1)   htemp=Xw_maineffects_self(x=xtemp[i], b=bvec, link=link, h.func=hfunc) else {
        htemp=Xw_maineffects_self(x=xtemp[i,], b=bvec, link=link, h.func=hfunc);
      };
      X.mat[i,]=htemp$X;
      w.vec[i]=htemp$w;
    };
    if(logscale) optemp=liftoneDoptimal_log_GLM_func (X=X.mat, w=w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram) else {
      optemp=liftoneDoptimal_GLM_func (X=X.mat, w=w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram);
    };
    m.design=sum(optemp$p>0);       # updated number of design point
    if(k.continuous==1) x.design=xtemp[optemp$p>0] else {
      x.design=xtemp[optemp$p>0,];  # updated list of design points
    };
    p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
    det.design=optemp$Maximum;      # optimized |X'WX|
    X.mat = X.mat[optemp$p>0,];   # updated model matrix
    w.vec = w.vec[optemp$p>0];    # updated w vector
    Dmat = t(X.mat * (p.design*w.vec)) %*% X.mat;    # X^T W X
    Amat = solve(Dmat);           # (X^T W X)^{-1}
    Qy <- function(y) {           # -Q(y1) given y1
      hy=hfunc(y);                # h(y)
      -(hy%*%Amat%*%hy)[1]*nutemp(sum(bvec*hy));
    };
    gradQ <- function(y) {        # gradient of -Q(y1)
      hy=hfunc(y);                # h(y)
      Ahy=(Amat%*%hy)[,1];        # gamma = A h(y)
      bhy=sum(bvec*hy);           # beta^T h(y)
      -nu1temp(bhy)*sum(hy*Ahy)*bvec[1:k.continuous] - 2*nutemp(bhy)*Ahy[1:k.continuous];
    };
    x0=(lvec+uvec)/2;
    ytemp=stats::optim(par=x0, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec,
                control=list(maxit=maxit, factr=reltol*1e13)); #(i)
    if(random) for(ia in 1:nram) { #random initial point, repeat (i), find potential better ytemp
      x0r=x0;
      for(i in 1:d.factor) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
      ytemp1=stats::optim(par=x0r, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec,
                   control=list(maxit=maxit, factr=reltol*1e13));
      if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
    };
    nit=1;                        # number of candidate y searched
    while((-ytemp$value/p.factor-1 > reltol)&&(nit < maxit)) { #start of step(ii)
      ystar=ytemp$par;                # candidate y
      hystar=hfunc(ystar);            # h(ystar)
      wstar=nutemp(sum(bvec*hystar)); # nu(beta^T h(y))
      alphat=0;                       # calculate weight of ystar
      bty=det.design;                                 # b_t
      dty=det(Dmat/2 + wstar/2*hystar%*%t(hystar));   # d_t
      py=p.factor;                                    # p
      if(2^py*dty > (py+1)*bty) {
        alphat=(2^py*dty - (py+1)*bty)/(py*(2^py*dty-2*bty));  # alpha_t
      };
      dtemp=function(x) { sqrt(sum((hystar-x)^2)); };
      atemp=apply(X.mat, 1, dtemp);
      #merge in(ii)
      if(min(atemp)<rel.diff) {   # merge "ystar" into its 1st neighbor
        iy=which.min(atemp);                          # index of design point to be merged
        p.design=(1-alphat)*p.design;
        if(k.continuous==1) {
          ystar1=(x.design[iy]*p.design[iy]+ystar*alphat)/(p.design[iy]+alphat);
          x.design[iy]=ystar1;                           # update x.design
        } else {
          ystar1=(x.design[iy,]*p.design[iy]+ystar*alphat)/(p.design[iy]+alphat);  # merged design point
          x.design[iy,]=ystar1;                           # update x.design
        };
        p.design[iy]=p.design[iy]+alphat;               # update p.design
        hystar1=hfunc(ystar1);                          # update X.mat
        X.mat[iy,]=hystar1;                             # update X.mat
        w.vec[iy]=nutemp(sum(bvec*hystar1));            # update w.vec
      } else {
        if(k.continuous==1) x.design=c(x.design,ystar) else {
          x.design=rbind(x.design,ystar);  # updated list of design points
        };
        X.mat=rbind(X.mat,hystar);    # add h(y) into design matrix
        w.vec=c(w.vec,wstar);
        p.design=c((1-alphat)*p.design, alphat);
      };                            # end of "if(...reltol)"
      if(logscale) optemp=liftoneDoptimal_log_GLM_func (X=X.mat, w=w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram) else {
        optemp=liftoneDoptimal_GLM_func (X=X.mat, w=w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram);
      };
      m.design=sum(optemp$p>0);       # updated number of design point
      if(k.continuous==1) x.design=x.design[optemp$p>0] else {
        x.design=x.design[optemp$p>0,];  # updated list of design points
      };
      p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
      det.design=optemp$Maximum;      # optimized |X'WX|
      X.mat = X.mat[optemp$p>0,];   # updated model matrix
      w.vec = w.vec[optemp$p>0];    # updated w vector
      Dmat = t(X.mat * (p.design*w.vec)) %*% X.mat;    # X^T W X
      Amat = solve(Dmat);           # (X^T W X)^{-1}
      Qy <- function(y) {           # -Q(y1) given y1
        hy=hfunc(y);                # h(y)
        -(hy%*%Amat%*%hy)[1]*nutemp(sum(bvec*hy));
      };
      gradQ <- function(y) {        # gradient of -Q(y1)
        hy=hfunc(y);                # h(y)
        Ahy=(Amat%*%hy)[,1];        # gamma = A h(y)
        bhy=sum(bvec*hy);           # beta^T h(y)
        -nu1temp(bhy)*sum(hy*Ahy)*bvec[1:k.continuous] - 2*nutemp(bhy)*Ahy[1:k.continuous];
      };
      x0=(lvec+uvec)/2;
      ytemp=stats::optim(par=x0, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec,
                  control=list(maxit=maxit, factr=reltol*1e13));
      if(random) for(ia in 1:nram) {
        x0r=x0;
        for(i in 1:d.factor) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
        ytemp1=stats::optim(par=x0r, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec,
                     control=list(maxit=maxit, factr=reltol*1e13));
        if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
      };
      nit=nit+1;                    # number of candidate y searched
    };     # End of "while" loop
    itmax.design=nit;
    converge.design=(ytemp$convergence==0);  # TRUE or FALSE
    if(-ytemp$value/p.factor-1 > reltol) converge.design=FALSE;
  };                                # end of Case II
  #    Case III: some factors are continuous
  if((k.continuous>0)&&(k.continuous<d.factor)) {
    lvec=uvec=rep(0, k.continuous);     # lower bounds and upper founds for continuous factors
    for(i in 1:k.continuous) {lvec[i]=min(factor.level[[i]]); uvec[i]=max(factor.level[[i]]);};
    if(is.null(Xini)) xtemp=xmat_discrete_self(factor.level, rowmax=rowmax) else xtemp=Xini;
    m.design=nrow(xtemp);                   # initial number of design points
    X.mat = matrix(0, m.design, p.factor);  # initial model matrix X
    w.vec = rep(0, m.design);     # w vector
    for(i in 1:m.design) {
      htemp=Xw_maineffects_self(x=xtemp[i,], b=bvec, link=link, h.func=hfunc);
      X.mat[i,]=htemp$X;
      w.vec[i]=htemp$w;
    };
    if(logscale) optemp=liftoneDoptimal_log_GLM_func (X=X.mat, w=w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram) else {
      optemp=liftoneDoptimal_GLM_func (X=X.mat, w=w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram);
    };
    m.design=sum(optemp$p>0);       # updated number of design point
    x.design=xtemp[optemp$p>0,]  # updated list of design points
    p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
    det.design=optemp$Maximum;      # optimized |X'WX|
    X.mat = X.mat[optemp$p>0,];   # updated model matrix
    w.vec = w.vec[optemp$p>0];    # updated w vector
    Dmat = t(X.mat * (p.design*w.vec)) %*% X.mat;    # X^T W X
    Amat = svd_inverse(Dmat);           # (X^T W X)^{-1}
    xdiscrete=xmat_discrete_self(factor.level[(k.continuous+1):d.factor]);
    ndiscrete=dim(xdiscrete)[1];
    for(idiscrete in 1:ndiscrete) {
      hfunc1 <- function(y) { hfunc(c(y, xdiscrete[idiscrete,])); };
      Qy <- function(y) {           # -Q(y1) given y1
        hy=hfunc1(y);               # h(y)
        -(hy%*%Amat%*%hy)[1]*nutemp(sum(bvec*hy));
      };
      gradQ <- function(y) {        # gradient of -Q(y1)
        hy=hfunc1(y);               # h(y)
        Ahy=(Amat%*%hy)[,1];        # gamma = A h(y)
        bhy=sum(bvec*hy);           # beta^T h(y)
        -nu1temp(bhy)*sum(hy*Ahy)*bvec[1:k.continuous] - 2*nutemp(bhy)*Ahy[1:k.continuous];
      };
      x0=(lvec+uvec)/2;
      ytemp=stats::optim(par=x0, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec,
                  control=list(maxit=maxit, factr=reltol*1e13));
      if(random) for(ia in 1:nram) {
        x0r=x0;
        for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
        ytemp1=stats::optim(par=x0r, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec,
                     control=list(maxit=maxit, factr=reltol*1e13));
        if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
      };
      if(idiscrete==1) {
        ystar=c(ytemp$par, xdiscrete[idiscrete,]);
        fvalue=ytemp$value;
      } else if(ytemp$value<fvalue) {
        ystar=c(ytemp$par, xdiscrete[idiscrete,]);
        fvalue=ytemp$value;
      };
    };                            # end of "idiscrete" loop
    nit=1;                        # number of candidate y searched
    while((-fvalue/p.factor-1 > reltol)&&(nit < maxit)) {
      hystar=hfunc(ystar);            # h(ystar)
      wstar=nutemp(sum(bvec*hystar)); # nu(beta^T h(y))
      alphat=0;                       # calculate weight of ystar
      bty=det.design;                                 # b_t
      dty=det(Dmat/2 + wstar/2*hystar%*%t(hystar));   # d_t
      py=p.factor;                                    # p
      if(2^py*dty > (py+1)*bty) {
        alphat=(2^py*dty - (py+1)*bty)/(py*(2^py*dty-2*bty));  # alpha_t
      };
      dtemp=function(x) { sqrt(sum((hystar-x)^2)); };
      atemp=apply(X.mat, 1, dtemp);
      if(min(atemp)<rel.diff) {   # merge "ystar" into its 1st neighbor
        iy=which.min(atemp);                          # index of design point to be merged
        p.design=(1-alphat)*p.design;
        ystar1=(x.design[iy,]*p.design[iy]+ystar*alphat)/(p.design[iy]+alphat);  # merged design point
        x.design[iy,]=ystar1;                           # update x.design
        p.design[iy]=p.design[iy]+alphat;               # update p.design
        hystar1=hfunc(ystar1);                         # update X.mat
        X.mat[iy,]=hystar1;                             # update X.mat
        w.vec[iy]=nutemp(sum(bvec*hystar1));            # update w.vec
      } else {
        x.design=rbind(x.design,ystar);  # updated list of design points
        X.mat=rbind(X.mat,hystar);    # add h(y) into design matrix
        w.vec=c(w.vec,wstar);
        p.design=c((1-alphat)*p.design, alphat);
      };                            # end of "if(...reltol)"
      if(logscale) optemp=liftoneDoptimal_log_GLM_func (X=X.mat, w=w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram) else {
        optemp=liftoneDoptimal_GLM_func (X=X.mat, w=w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram);
      };
      m.design=sum(optemp$p>0);       # updated number of design point
      x.design=x.design[optemp$p>0,];  # updated list of design points
      p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
      det.design=optemp$Maximum;      # optimized |X'WX|
      X.mat = X.mat[optemp$p>0,];   # updated model matrix
      w.vec = w.vec[optemp$p>0];    # updated w vector
      Dmat = t(X.mat * (p.design*w.vec)) %*% X.mat;    # X^T W X
      Amat = svd_inverse(Dmat);           # (X^T W X)^{-1}
      for(idiscrete in 1:ndiscrete) {
        hfunc1 <- function(y) { hfunc(c(y, xdiscrete[idiscrete,])); };
        Qy <- function(y) {           # -Q(y1) given y1
          hy=hfunc1(y);               # h(y)
          -(hy%*%Amat%*%hy)[1]*nutemp(sum(bvec*hy));
        };
        gradQ <- function(y) {        # gradient of -Q(y1)
          hy=hfunc1(y);               # h(y)
          Ahy=(Amat%*%hy)[,1];        # gamma = A h(y)
          bhy=sum(bvec*hy);           # beta^T h(y)
          -nu1temp(bhy)*sum(hy*Ahy)*bvec[1:k.continuous] - 2*nutemp(bhy)*Ahy[1:k.continuous];
        };
        x0=(lvec+uvec)/2;
        ytemp=stats::optim(par=x0, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec,
                    control=list(maxit=maxit, factr=reltol*1e13));
        if(random) for(ia in 1:nram) {
          x0r=x0;
          for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
          ytemp1=stats::optim(par=x0r, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec,
                       control=list(maxit=maxit, factr=reltol*1e13));
          if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
        };
        if(idiscrete==1) {
          ystar=c(ytemp$par, xdiscrete[idiscrete,]);
          fvalue=ytemp$value;
        } else if(ytemp$value<fvalue) {
          ystar=c(ytemp$par, xdiscrete[idiscrete,]);
          fvalue=ytemp$value;
        };
      };                            # end of "idiscrete" loop
      nit=nit+1;                    # number of candidate y searched
    };     # End of "while" loop
    itmax.design=nit;
    converge.design=(ytemp$convergence==0);  # TRUE or FALSE
    if(-ytemp$value/p.factor-1 > reltol) converge.design=FALSE;
  };                                # end of Case III
  rownames(x.design)=NULL;
  rownames(X.mat)=NULL;
  dtemp=as.matrix(stats::dist(x.design));
  diag(dtemp)=Inf;
  min.diff=min(dtemp);
  i1=which.min(apply(dtemp,1,min));
  i2=which.min(dtemp[i1,]);
  if(d.factor==1) x.close=x.design[c(i1,i2)] else {
    x.close=x.design[c(i1,i2),];
  };
  # list(m=m.design, x.factor=x.design, p=p.design, det=det.design, x.model=X.mat, w=w.vec, convergence=converge.design, min.diff=min.diff, x.close=x.close);

  #define S3 class
  output<-list(m=m.design, x.factor=x.design, p=p.design, det=det.design, convergence=converge.design, min.diff=min.diff, x.close=x.close, itmax=itmax.design);
  class(output) <- "design_output"
  return(output)
}









