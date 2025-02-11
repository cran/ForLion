#' EW ForLion for generalized linear models
#' @description
#' EW ForLion algorithm to find EW D-optimal design for GLM models with mixed factors, reference: .
#' Factors may include discrete factors with finite number of distinct levels and continuous factors
#' with specified interval range (min, max), continuous factors, if any, must serve as main-effects
#' only, allowing merging points that are close enough.Continuous factors first then discrete factors,
#'  model parameters should in the same order of factors.
#' @param n.factor vector of numbers of distinct levels, "0" indicates continuous factors, "0"s always come first, "2" or above indicates discrete factor, "1" is not allowed
#' @param factor.level list of distinct levels, (min, max) for continuous factor, continuous factors first, should be the same order as n.factor
#' @param hfunc function for obtaining model matrix h(y) for given design point y, y has to follow the same order as n.factor
#' @param joint_Func_b The prior joint probability distribution of the parameters
#' @param Lowerbounds The lower limit of the prior distribution for each parameter
#' @param Upperbounds The upper limit of the prior distribution for each parameter
#' @param link link function, default "logit", other links: "probit", "cloglog", "loglog", "cauchit", "log", "identity"
#' @param reltol the relative convergence tolerance, default value 1e-5
#' @param rel.diff points with distance less than that will be merged, default value 0
#' @param optim_grad TRUE or FALSE, default is FALSE, whether to use the analytical gradient function or numerical gradient for searching optimal new design point
#' @param maxit the maximum number of iterations, default value 100
#' @param random TRUE or FALSE, if TRUE then the function will run EW lift-one with additional "nram" number of random approximate allocation, default to be FALSE
#' @param nram when random == TRUE, the function will run EW lift-one nram number of initial proportion p00, default is 3
#' @param logscale TRUE or FALSE, if TRUE then the EW ForLion will run EW lift-one with logscale, which is EW_liftoneDoptimal_log_GLM_func(); if FALSE then ForLion will run EW lift-one without logscale, which is EW_liftoneDoptimal_GLM_func()
#' @param rowmax maximum number of points in the initial design, default NULL indicates no restriction
#' @param Xini initial list of design points, default NULL will generate random initial design points
#'
#' @return m number of design points
#' @return x.factor matrix with rows indicating design point
#' @return p EW D-optimal approximate allocation
#' @return det Optimal determinant of Fisher information matrix
#' @return x.model model matrix X
#' @return E_w vector of E_w such that E_w=diag(p*E_w)
#' @return convergence TRUE or FALSE
#' @return min.diff the minimum Euclidean distance between design points
#' @return x.close  a pair of design points with minimum distance
#' @export
#'
#' @examples
#' #Example  Crystallography Experiment
#' hfunc.temp = function(y) {c(y,1)}   # y -> h(y)=(y1,1)
#' n.factor.temp = c(0)  # 1 continuous factors
#' factor.level.temp = list(c(-1,1))
#' link.temp="logit"
#' paras_lowerbound<-c(4,-3)
#' paras_upperbound<-c(10,3)
#'  gjoint_b<- function(x) {
#'  Func_b<-1/(prod(paras_upperbound-paras_lowerbound))
#'  ##the prior distributions are follow uniform distribution
#' return(Func_b)
#' }
#' EW_ForLion_GLM_Optimal(n.factor=n.factor.temp, factor.level=factor.level.temp,
#' hfunc=hfunc.temp,joint_Func_b=gjoint_b, Lowerbounds=paras_lowerbound,
#' Upperbounds=paras_upperbound, link=link.temp, reltol=1e-2, rel.diff=0.01,
#' optim_grad=FALSE, maxit=500, random=FALSE, nram=3, logscale=FALSE,Xini=NULL)

EW_ForLion_GLM_Optimal<- function(n.factor, factor.level, hfunc,joint_Func_b,
          Lowerbounds, Upperbounds, link, reltol=1e-5, rel.diff=0,optim_grad=TRUE,
          maxit=100, random=FALSE, nram=3, logscale=FALSE, rowmax=NULL, Xini=NULL) {
  d.factor=length(n.factor);             # number of factors
  p.factor=length(Lowerbounds);                 # number of predictors
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
    E_w.vec = rep(0, m.design);     # E_w vector
    for(i in 1:m.design) {
      htemp=EW_Xw_maineffects_self(x=xtemp[i,],joint_Func_b, Lowerbounds, Upperbounds, link=link, h.func=hfunc);
      X.mat[i,]=htemp$X;
      E_w.vec[i]=htemp$E_w;
    };
    if(logscale) optemp=EW_liftoneDoptimal_log_GLM_func(X=X.mat, E_w=E_w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram) else {
      optemp=EW_liftoneDoptimal_GLM_func(X=X.mat, E_w=E_w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram);
    };
    m.design=sum(optemp$p>0);       # updated number of design point
    x.design=xtemp[optemp$p>0,];  # updated list of design points
    p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
    det.design=optemp$Maximum;      # optimized |X'E_wX|
    X.mat = X.mat[optemp$p>0,];   # updated model matrix
    E_w.vec = E_w.vec[optemp$p>0];    # updated E_w vector
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
    E_w.vec = rep(0, m.design);     # E_w vector
    for(i in 1:m.design) {
      if(k.continuous==1)   htemp=EW_Xw_maineffects_self(x=xtemp[i],joint_Func_b,Lowerbounds, Upperbounds,link=link, h.func=hfunc) else {
        htemp=EW_Xw_maineffects_self(x=xtemp[i,],joint_Func_b, Lowerbounds, Upperbounds, link=link, h.func=hfunc);
      };
      X.mat[i,]=htemp$X;
      E_w.vec[i]=htemp$E_w;
    };
    if(logscale) optemp=EW_liftoneDoptimal_log_GLM_func(X=X.mat, E_w=E_w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram) else {
      optemp=EW_liftoneDoptimal_GLM_func(X=X.mat, E_w=E_w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram);
    };
    m.design=sum(optemp$p>0);       # updated number of design point
    if(k.continuous==1) x.design=xtemp[optemp$p>0] else {
      x.design=xtemp[optemp$p>0,];  # updated list of design points
    };
    p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
    det.design=optemp$Maximum;      # optimized |X'E_wX|
    X.mat = X.mat[optemp$p>0,];   # updated model matrix
    E_w.vec = E_w.vec[optemp$p>0];    # updated E_w vector
    Dmat = t(X.mat * (p.design*E_w.vec)) %*% X.mat;    # X^T E_w X
    Amat = solve(Dmat);           # (X^T E_w X)^{-1}
    Qy <- function(y) {           # -Q(y1) given y1
      hy=hfunc(y); # h(y)
      integrand_dsensi<-function(b){
        etay = sum(b*hy);
        w_hy<-nutemp(etay)
        return(w_hy * joint_Func_b(b))
      }
      result_hy <- cubature::hcubature(f =integrand_dsensi,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
      Ewy<-result_hy$integral
      Ew_d<-(-(hy%*%Amat%*%hy)[1]*Ewy);
      return(Ew_d)
    };
    gradQ <- function(y) {# gradient of -Q(y1)
      Ew_d_dx<-rep(0,length(k.continuous))
      for(hq in 1:k.continuous){
        hy=hfunc(y);                # h(y)
        Ahy=(Amat%*%hy)[,1];# gamma = A h(y)
        integrand_d_dx<-function(b){
          bhy=sum(b*hy);           # beta^T h(y)
          d_dx_inter<- (-nu1temp(bhy)*sum(hy*Ahy)*b[hq] - 2*nutemp(bhy)*Ahy[hq]);
          return(d_dx_inter * joint_Func_b(b))
        }
        result_qd_dx <- cubature::hcubature(f =integrand_d_dx,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
        Ew_d_dx[hq]<-result_qd_dx$integral
      }
      return(Ew_d_dx)
    };
    x0=(lvec+uvec)/2;
    if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=Qy, gr=gradQ,method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
    if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=Qy, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}

    if(random) for(ia in 1:nram) { #random initial point, repeat (i), find potential better ytemp
      x0r=x0;
      for(i in 1:d.factor) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
      if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
      if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=Qy, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}

      if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
    };
    nit=1;                        # number of candidate y searched
    while((-ytemp$value/p.factor-1 > reltol)&&(nit < maxit)) { #start of step(ii)
      ystar=ytemp$par;                # candidate y
      hystar=hfunc(ystar);            # h(ystar)

      integrand_wstar<-function(b){
        etay_star = sum(b*hystar);
        w_hystar<-nutemp(etay_star)
        return(w_hystar * joint_Func_b(b))
      }
      result_wstar <- cubature::hcubature(f =integrand_wstar,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
      wstar<-result_wstar$integral  # nu(beta^T h(y))

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

        integrand_wstar1<-function(b){
          etay_star1 = sum(b*hystar1);
          w_hystar1<-nutemp(etay_star1)
          return(w_hystar1 * joint_Func_b(b))
        }
        result_wstar1 <- cubature::hcubature(f =integrand_wstar1,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
        E_w.vec[iy]=result_wstar1$integral   # update w.vec
      } else {
        if(k.continuous==1) x.design=c(x.design,ystar) else {
          x.design=rbind(x.design,ystar);  # updated list of design points
        };
        X.mat=rbind(X.mat,hystar);    # add h(y) into design matrix
        E_w.vec=c(E_w.vec,wstar);
        p.design=c((1-alphat)*p.design, alphat);
      };                            # end of "if(...reltol)"
      if(logscale) optemp=EW_liftoneDoptimal_log_GLM_func(X=X.mat, E_w=E_w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram) else {
        optemp=EW_liftoneDoptimal_GLM_func(X=X.mat, E_w=E_w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram);
      };
      m.design=sum(optemp$p>0);       # updated number of design point
      if(k.continuous==1) x.design=x.design[optemp$p>0] else {
        x.design=x.design[optemp$p>0,];  # updated list of design points
      };
      p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
      det.design=optemp$Maximum;      # optimized |X'E_wX|
      X.mat = X.mat[optemp$p>0,];   # updated model matrix
      E_w.vec = E_w.vec[optemp$p>0];    # updated E_w vector
      Dmat = t(X.mat * (p.design*E_w.vec)) %*% X.mat;    # X^T E_w X
      Amat = solve(Dmat);           # (X^T E_w X)^{-1}
      Qy <- function(y) {           # -Q(y1) given y1
        hy=hfunc(y); # h(y)
        integrand_dsensi<-function(b){
          etay = sum(b*hy);
          w_hy<-nutemp(etay)
          return(w_hy * joint_Func_b(b))
        }
        result_hy <- cubature::hcubature(f =integrand_dsensi,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
        Ewy<-result_hy$integral
        Ew_d<-(-(hy%*%Amat%*%hy)[1]*Ewy);
        return(Ew_d)
      };
      gradQ <- function(y) {# gradient of -Q(y1)
        Ew_d_dx<-NULL
        for(hq in 1:k.continuous){
          hy=hfunc(y);                # h(y)
          Ahy=(Amat%*%hy)[,1];# gamma = A h(y)
          integrand_d_dx<-function(b){
            bhy=sum(b*hy);           # beta^T h(y)
            d_dx_inter<- (-nu1temp(bhy)*sum(hy*Ahy)*b[hq] - 2*nutemp(bhy)*Ahy[hq]);
            return(d_dx_inter * joint_Func_b(b))
          }
          result_qd_dx <- cubature::hcubature(f =integrand_d_dx,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
          Ew_d_dx[hq]<-result_qd_dx$integral
        }
        return(Ew_d_dx)
      };

      x0=(lvec+uvec)/2;
      if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
      if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=Qy, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}

      if(random) for(ia in 1:nram) {
        x0r=x0;
        for(i in 1:d.factor) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
        if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
        if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=Qy, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}

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
    E_w.vec = rep(0, m.design);     # E_w vector
    for(i in 1:m.design) {
      htemp=EW_Xw_maineffects_self(x=xtemp[i, ],joint_Func_b,Lowerbounds, Upperbounds, link=link, h.func=hfunc);
      X.mat[i,]=htemp$X;
      E_w.vec[i]=htemp$E_w;
    };
    if(logscale) optemp=EW_liftoneDoptimal_log_GLM_func(X=X.mat, E_w=E_w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram) else {
      optemp=EW_liftoneDoptimal_GLM_func(X=X.mat, E_w=E_w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram);
    };
    m.design=sum(optemp$p>0);       # updated number of design point
    x.design=xtemp[optemp$p>0,]  # updated list of design points
    p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
    det.design=optemp$Maximum;      # optimized |X'E_wX|
    X.mat = X.mat[optemp$p>0,];   # updated model matrix
    E_w.vec = E_w.vec[optemp$p>0];    # updated E_w vector
    Dmat = t(X.mat * (p.design*E_w.vec)) %*% X.mat;    # X^T E_w X
    Amat = svd_inverse(Dmat);           # (X^T E_w X)^{-1}
    xdiscrete=xmat_discrete_self(factor.level[(k.continuous+1):d.factor]);
    ndiscrete=dim(xdiscrete)[1];
    for(idiscrete in 1:ndiscrete) {
      hfunc1 <- function(y) { hfunc(c(y, xdiscrete[idiscrete,])); };
      Qy <- function(y) {           # -Q(y1) given y1
        hy=hfunc1(y); # h(y)
        integrand_dsensi<-function(b){
          etay = sum(b*hy);
          w_hy<-nutemp(etay)
          return(w_hy * joint_Func_b(b))
        }
        result_hy <- cubature::hcubature(f =integrand_dsensi,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
        Ewy<-result_hy$integral
        Ew_d<-(-(hy%*%Amat%*%hy)[1]*Ewy);
        return(Ew_d)
      };
      gradQ <- function(y) {# gradient of -Q(y1)
        Ew_d_dx<-rep(0,length(k.continuous))
        for(hq in 1:k.continuous){
          hy=hfunc1(y);                # h(y)
          Ahy=(Amat%*%hy)[,1];# gamma = A h(y)
          integrand_d_dx<-function(b){
            bhy=sum(b*hy);           # beta^T h(y)
            d_dx_inter<- (-nu1temp(bhy)*sum(hy*Ahy)*b[hq] - 2*nutemp(bhy)*Ahy[hq]);
            return(d_dx_inter * joint_Func_b(b))
          }
          result_qd_dx <- cubature::hcubature(f =integrand_d_dx,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
          Ew_d_dx[hq]<-result_qd_dx$integral
        }
        return(Ew_d_dx)
      };
      x0=(lvec+uvec)/2;
      if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=Qy, gr=gradQ,method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
      if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=Qy, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}

      if(random) for(ia in 1:nram) {
        x0r=x0;
        for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
        if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
        if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=Qy, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}

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

      integrand_wstar<-function(b){
        etay_star = sum(b*hystar);
        w_hystar<-nutemp(etay_star)
        return(w_hystar * joint_Func_b(b))
      }
      result_wstar <- cubature::hcubature(f =integrand_wstar,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
      wstar<-result_wstar$integral  # nu(beta^T h(y))

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

        integrand_wstar1<-function(b){
          etay_star1 = sum(b*hystar1);
          w_hystar1<-nutemp(etay_star1)
          return(w_hystar1 * joint_Func_b(b))
        }
        result_wstar1 <- cubature::hcubature(f =integrand_wstar1,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
        E_w.vec[iy]=result_wstar1$integral   # update w.vec
      } else {
        x.design=rbind(x.design,ystar);  # updated list of design points
        X.mat=rbind(X.mat,hystar);    # add h(y) into design matrix
        E_w.vec=c(E_w.vec,wstar);
        p.design=c((1-alphat)*p.design, alphat);
      };                            # end of "if(...reltol)"
      if(logscale) optemp=EW_liftoneDoptimal_log_GLM_func(X=X.mat, E_w=E_w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram) else {
        optemp=EW_liftoneDoptimal_GLM_func(X=X.mat, E_w=E_w.vec, reltol=reltol, maxit=maxit, random=random, nram=nram);
      };
      m.design=sum(optemp$p>0);       # updated number of design point
      x.design=x.design[optemp$p>0,];  # updated list of design points
      p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
      det.design=optemp$Maximum;      # optimized |X'E_wX|
      X.mat = X.mat[optemp$p>0,];   # updated model matrix
      E_w.vec = E_w.vec[optemp$p>0];    # updated E_w vector
      Dmat = t(X.mat * (p.design*E_w.vec)) %*% X.mat;    # X^T E_w X
      Amat = svd_inverse(Dmat);           # (X^T E_w X)^{-1}
      for(idiscrete in 1:ndiscrete) {
        hfunc1 <- function(y) { hfunc(c(y, xdiscrete[idiscrete,])); };
        Qy <- function(y) {           # -Q(y1) given y1
          hy=hfunc1(y); # h(y)
          integrand_dsensi<-function(b){
            etay = sum(b*hy);
            w_hy<-nutemp(etay)
            return(w_hy * joint_Func_b(b))
          }
          result_hy <- cubature::hcubature(f =integrand_dsensi,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
          Ewy<-result_hy$integral
          Ew_d<-(-(hy%*%Amat%*%hy)[1]*Ewy);
          return(Ew_d)
        };
        gradQ <- function(y) {# gradient of -Q(y1)
          Ew_d_dx<-rep(0,length(k.continuous))
          for(hq in 1:k.continuous){
            hy=hfunc1(y);                # h(y)
            Ahy=(Amat%*%hy)[,1];# gamma = A h(y)
            integrand_d_dx<-function(b){
              bhy=sum(b*hy);           # beta^T h(y)
              d_dx_inter<- (-nu1temp(bhy)*sum(hy*Ahy)*b[hq] - 2*nutemp(bhy)*Ahy[hq]);
              return(d_dx_inter * joint_Func_b(b))
            }
            result_qd_dx <- cubature::hcubature(f =integrand_d_dx,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
            Ew_d_dx[hq]<-result_qd_dx$integral
          }
          return(Ew_d_dx)
        };
        x0=(lvec+uvec)/2;
        if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=Qy, gr=gradQ,method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
        if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=Qy, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}

        if(random) for(ia in 1:nram) {
          x0r=x0;
          for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
          if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
          if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=Qy, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}

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
  #list(m=m.design, x.factor=x.design, p=p.design, det=det.design, x.model=X.mat,
  #     E_w=E_w.vec, convergence=converge.design, min.diff=min.diff, x.close=x.close);

  #define S3 class
  output<-list(m=m.design, x.factor=x.design, p=p.design, det=det.design, convergence=converge.design, min.diff=min.diff, x.close=x.close, itmax=itmax.design);
  class(output) <- "design_output"
  return(output)
}



