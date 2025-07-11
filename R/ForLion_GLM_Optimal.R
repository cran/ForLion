#' ForLion for generalized linear models
#' @description
#' ForLion algorithm to find D-optimal design for GLM models with mixed factors, reference: Section 4 in Huang, Li, Mandal, Yang (2024).
#' Factors may include discrete factors with finite number of distinct levels and continuous factors with specified interval range (min, max), continuous factors, if any, must serve as main-effects only, allowing merging points that are close enough.
#' Continuous factors first then discrete factors, model parameters should in the same order of factors.
#' @param n.factor vector of numbers of distinct levels, “0” indicating continuous factors that always come first, “2” or more for discrete factors, and “1” not allowed.
#' @param factor.level list of distinct factor levels, “(min, max)” for continuous factors that always come first, finite sets for discrete factors.
#' @param var_names Names for the design factors. Must have the same length asfactor.level. Defaults to "X1", "X2", ...
#' @param xlist_fix list of discrete factor experimental settings under consideration, default NULL indicating a list of all possible discrete factor experimental settings will be used.
#' @param hfunc function for generating the corresponding model matrix or predictor vector, given an experimental setting or design point.
#' @param bvec assumed parameter values of model parameters beta, same length of h(y)
#' @param link link function, default "logit", other links: "probit", "cloglog", "loglog", "cauchit", "log", "identity"
#' @param reltol the relative convergence tolerance, default value 1e-5
#' @param delta relative difference as merging threshold for the merging step, the distance of two points less than delta may be merged, default 0, can be different from delta0 for the initial design.
#' @param maxit the maximum number of iterations, default value 100
#' @param random TRUE or FALSE, if TRUE then the function will run lift-one with additional "nram" number of random approximate allocation, default to be FALSE
#' @param nram when random == TRUE, the function will run lift-one nram number of initial proportion p00, default is 3
#' @param logscale TRUE or FALSE, whether or not to run the lift-one step in log-scale, i.e., using liftoneDoptimal_log_GLM_func() or liftoneDoptimal_GLM_func().
#' @param rowmax maximum number of points in the initial design, default NULL indicates no restriction
#' @param Xini initial list of design points, default NULL indicating automatically generating an initial list of design points.
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
#' ForLion_GLM_Optimal(n.factor=n.factor.temp, factor.level=factor.level.temp, xlist_fix=NULL,
#' hfunc=hfunc.temp, bvec=b.temp, link=link.temp, reltol=1e-2, delta=0.03, maxit=500,
#' random=FALSE,nram=3, logscale=TRUE)
#'


ForLion_GLM_Optimal <- function(n.factor, factor.level,var_names=NULL,xlist_fix=NULL, hfunc, bvec, link, reltol=1e-5, delta=0, maxit=100, random=FALSE, nram=3, logscale=FALSE, rowmax=NULL, Xini=NULL) {
  d.factor=length(n.factor);             # number of factors
  p.factor=length(bvec);                 # number of predictors
  k.continuous=sum(n.factor==0);         # number of continuous factors
  if(delta==0) delta=reltol;
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
    if(is.null(Xini)){initial.temp=design_initial_self(k.continuous=k.continuous, factor.level=factor.level, MLM=FALSE, xlist_fix=xlist_fix, lvec=lvec, uvec=uvec, bvec=bvec, link=link, h.func=hfunc, delta0=reltol, epsilon = reltol, maxit=maxit); xtemp=initial.temp$X} else {xtemp=Xini}  #no initial design
    if(k.continuous==1) m.design=length(xtemp) else m.design=nrow(xtemp);                   # initial number of design points
    X.mat = matrix(0, m.design, p.factor);  # initial model matrix X
    w.vec = rep(0, m.design);               # w vector
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

      #add new points to the design
      if(k.continuous==1) x.design=c(x.design,ystar) else {
        x.design=rbind(x.design,ystar);  # updated list of design points
      };
      X.mat=rbind(X.mat,hystar);    # add h(y) into design matrix
      w.vec=c(w.vec,wstar);
      p.design=c((1-alphat)*p.design, alphat);

      #merging: calculate distance between design points, if min distance < tolerance merge
      dtemp=as.matrix(stats::dist(x.design));
      diag(dtemp)=Inf;
      atemp=min(dtemp)
      while((atemp<delta)){ # merge closest two neighbors
        #before merging two closest points, save the current state of the design
        x.design_old=x.design
        p.design_old=p.design
        X.mat_old=X.mat
        w.vec_old=w.vec

        #identify and merge the two closest design points
        i1=which.min(apply(dtemp,1,min)); # index of design point to be merged
        i2=which.min(dtemp[i1,]); # index of design point to be merged
        if(k.continuous==1){
          ystar1=(p.design_old[i1]*x.design_old[i1]+p.design_old[i2]*x.design_old[i2])/(p.design_old[i1]+p.design_old[i2]); # merged design point
          x.design_mer=c(x.design_old[-c(i1,i2)], ystar1)
        }else{
          ystar1=(p.design_old[i1]*x.design_old[i1,]+p.design_old[i2]*x.design_old[i2,])/(p.design_old[i1]+p.design_old[i2]); # merged design point
          x.design_mer=rbind(x.design_old[-c(i1,i2), ], ystar1) # update x.design
        }
        hystar1=hfunc(ystar1);            # h(ystar1)
        wstar1=nutemp(sum(bvec*hystar1)); # nu(beta^T h(y))

        X.mat_mer=rbind(X.mat_old[-c(i1,i2), ],hystar1);    # add h(y) into design matrix
        w.vec_mer=c(w.vec_old[-c(i1,i2)],wstar1);
        p.design_mer=c(p.design_old[-c(i1,i2)], p.design_old[i1]+p.design_old[i2])

        # Build X.mat_mer and calculate the Fisher information matrix (F.mat_mer)
        F.mat_mer=det(t(X.mat_mer * (p.design_mer*w.vec_mer)) %*% X.mat_mer)

        eigen_values<-eigen(F.mat_mer)
        min_engenvalue<-min(eigen_values$values)
        if(min_engenvalue<=reltol){
          x.design=x.design_old
          p.design=p.design_old
          X.mat=X.mat_old
          w.vec=w.vec_old
          break
        }else{
          x.design=x.design_mer
          p.design=p.design_mer
          X.mat=X.mat_mer
          w.vec=w.vec_mer
        }
        dtemp=as.matrix(stats::dist(x.design));
        diag(dtemp)=Inf;
        atemp=min(dtemp)
      }

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
    if(is.null(Xini)){initial.temp=design_initial_self(k.continuous=k.continuous, factor.level=factor.level, MLM=FALSE, xlist_fix=xlist_fix, lvec=lvec, uvec=uvec, bvec=bvec, link=link, h.func=hfunc, delta0=reltol, epsilon = reltol, maxit=maxit); xtemp=initial.temp$X} else {xtemp=Xini}  #no initial design
    #if(is.null(Xini)) xtemp=xmat_discrete_self(factor.level, rowmax=rowmax) else xtemp=Xini;
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
    if(is.null(xlist_fix)){xdiscrete=xmat_discrete_self(factor.level[(k.continuous+1):d.factor]);}else{xdiscrete=xlist_fix;}
    # xdiscrete=xmat_discrete_self(factor.level[(k.continuous+1):d.factor]);
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
      ytempstar=ytemp;
      if(random) for(ia in 1:nram) {
        x0r=x0;
        for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
        ytemp1=stats::optim(par=x0r, fn=Qy, gr=gradQ, method="L-BFGS-B", lower=lvec, upper=uvec,
                            control=list(maxit=maxit, factr=reltol*1e13));
        if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1;};
      };
      if(idiscrete==1) {
        ystar=c(ytemp$par, xdiscrete[idiscrete,]);
        fvalue=ytemp$value;
        ytempstar=ytemp;
      } else if(ytemp$value<fvalue) {
        ystar=c(ytemp$par, xdiscrete[idiscrete,]);
        fvalue=ytemp$value;
        ytempstar=ytemp;
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

      #add new points to the design
      x.design=rbind(x.design,ystar);  # updated list of design points
      X.mat=rbind(X.mat,hystar);    # add h(y) into design matrix
      w.vec=c(w.vec,wstar);
      p.design=c((1-alphat)*p.design, alphat);

      #merging: calculate distance between design points, if min distance < tolerance merge
      dtemp=as.matrix(stats::dist(x.design));
      diag(dtemp)=Inf;
      atemp=min(dtemp)
      while((atemp<delta)){ # merge closest two neighbors
        #before merging two closest points, save the current state of the design
        x.design_old=x.design
        p.design_old=p.design
        X.mat_old=X.mat
        w.vec_old=w.vec

        #identify and merge the two closest design points
        i1=which.min(apply(dtemp,1,min)); # index of design point to be merged
        i2=which.min(dtemp[i1,]); # index of design point to be merged

        ystar1=(p.design_old[i1]*x.design_old[i1,]+p.design_old[i2]*x.design_old[i2,])/(p.design_old[i1]+p.design_old[i2]); # merged design point
        x.design_mer=rbind(x.design_old[-c(i1,i2), ], ystar1) # update x.design

        hystar1=hfunc(ystar1);            # h(ystar1)
        wstar1=nutemp(sum(bvec*hystar1)); # nu(beta^T h(y))

        X.mat_mer=rbind(X.mat_old[-c(i1,i2), ],hystar1);    # add h(y) into design matrix
        w.vec_mer=c(w.vec_old[-c(i1,i2)],wstar1);
        p.design_mer=c(p.design_old[-c(i1,i2)], p.design_old[i1]+p.design_old[i2])

        # Build X.mat_mer and calculate the Fisher information matrix (F.mat_mer)
        F.mat_mer=det(t(X.mat_mer * (p.design_mer*w.vec_mer)) %*% X.mat_mer)

        eigen_values<-eigen(F.mat_mer)
        min_engenvalue<-min(eigen_values$values)

        if(min_engenvalue<=reltol){
          x.design=x.design_old
          p.design=p.design_old
          X.mat=X.mat_old
          w.vec=w.vec_old
          break
        }else{
          x.design=x.design_mer
          p.design=p.design_mer
          X.mat=X.mat_mer
          w.vec=w.vec_mer
        }
        dtemp=as.matrix(stats::dist(x.design));
        diag(dtemp)=Inf;
        atemp=min(dtemp)
      }

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
          if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1;};
        };
        if(idiscrete==1) {
          ystar=c(ytemp$par, xdiscrete[idiscrete,]);
          fvalue=ytemp$value;
          ytempstar=ytemp;
        } else if(ytemp$value<fvalue) {
          ystar=c(ytemp$par, xdiscrete[idiscrete,]);
          fvalue=ytemp$value;
          ytempstar=ytemp;
        };
      };                            # end of "idiscrete" loop
      nit=nit+1;                    # number of candidate y searched
    };     # End of "while" loop
    itmax.design=nit;
    converge.design=(ytempstar$convergence==0);  # TRUE or FALSE
    if(-ytempstar$value/p.factor-1 > reltol) converge.design=FALSE;
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
  output<-list(m=m.design, x.factor=x.design, p=p.design,var.names=var_names, det=det.design, convergence=converge.design, min.diff=min.diff, x.close=x.close, itmax=itmax.design);
  class(output) <- "design_output"
  return(output)
}
