#' ForLion function for multinomial logit models
#' @description
#' Function for ForLion algorithm to find D-optimal design under multinomial logit models with mixed factors.
#' Reference Section 3 of Huang, Li, Mandal, Yang (2024).
#' Factors may include discrete factors with finite number of distinct levels and continuous factors with specified interval range (min, max), continuous factors, if any, must serve as main-effects only, allowing merging points that are close enough.
#' Continuous factors first then discrete factors, model parameters should in the same order of factors.
#' @param J number of response levels in the multinomial logit model
#' @param n.factor Vector of numbers of distinct levels, “0” indicating continuous factors that always come first, “2” or more for discrete factors, and “1” not allowed.
#' @param factor.level list of distinct factor levels, “(min, max)” for continuous factors that always come first, finite sets for discrete factors.
#' @param var_names Names for the design factors. Must have the same length asfactor.level. Defaults to "X1", "X2", ...
#' @param xlist_fix list of discrete factor experimental settings under consideration, default NULL indicating a list of all possible discrete factor experimental settings will be used.
#' @param hfunc function for generating the corresponding model matrix or predictor vector, given an experimental setting or design point.
#' @param h.prime function to obtain dX/dx
#' @param bvec assumed parameter values of model parameters beta, same length of h(y)
#' @param link link function, default "continuation", other choices "baseline", "cumulative", and "adjacent"
#' @param Fi.func function to calculate row-wise Fisher information Fi, default is Fi_MLM_func
#' @param delta0  merging threshold for initial design, such that, || x_i(0) - x_j(0) || >= delta0, default 1e-5
#' @param epsilon tuning parameter as converging threshold, such that, a nonnegative number is regarded as numerical zero if less than epsilon, default 1e-12.
#' @param reltol the relative convergence tolerance, default value 1e-5
#' @param delta relative difference as merging threshold for the merging step, the distance of two points less than delta may be merged, default 0, can be different from delta0 for the initial design.
#' @param maxit the maximum number of iterations, default value 100
#' @param random TRUE or FALSE, whether or not to repeat the lift-one step multiple times with random initial allocations, default FALSE.
#' @param nram number of times repeating the lift-one step with random initial allocations, valid only if random is TRUE, default 3.
#' @param rowmax maximum number of points in the initial design, default NULL indicates no restriction
#' @param Xini initial list of design points, default NULL indicating automatically generating an initial list of design points.
#' @param random.initial TRUE or FALSE, whether or not to repeat the whole procedure multiple times with random initial designs, default FALSE.
#' @param nram.initial number of times repeating the whole procedure with random initial designs, valid only if random.initial is TRUE, default 3.
#' @param optim_grad TRUE or FALSE, default is FALSE, whether to use the analytical gradient function or numerical gradient when searching for a new design point.
#'
#' @return m             the number of design points
#' @return x.factor      matrix of experimental factors with rows indicating design point
#' @return p             the reported D-optimal approximate allocation
#' @return det           the determinant of the maximum Fisher information
#' @return convergence   TRUE or FALSE, whether converge
#' @return min.diff      the minimum Euclidean distance between design points
#' @return x.close       pair of design points with minimum distance
#' @return itmax         iteration of the algorithm
#' @export
#'
#' @examples
#' m=5
#' p=10
#' J=5
#' link.temp = "cumulative"
#' n.factor.temp = c(0,0,0,0,0,2)  # 1 discrete factor w/ 2 levels + 5 continuous
#' ## Note: Always put continuous factors ahead of discrete factors,
#' ## pay attention to the order of coefficients paring with predictors
#' factor.level.temp = list(c(-25,25), c(-200,200),c(-150,0),c(-100,0),c(0,16),c(-1,1))
#' hfunc.temp = function(y){
#' if(length(y) != 6){stop("Input should have length 6");}
#'  model.mat = matrix(NA, nrow=5, ncol=10, byrow=TRUE)
#'  model.mat[5,]=0
#'  model.mat[1:4,1:4] = diag(4)
#'  model.mat[1:4, 5] =((-1)*y[6])
#'  model.mat[1:4, 6:10] = matrix(((-1)*y[1:5]), nrow=4, ncol=5, byrow=TRUE)
#'  return(model.mat)
#'  }
#' bvec.temp=c(-1.77994301, -0.05287782,  1.86852211, 2.76330779, -0.94437464, 0.18504420,
#' -0.01638597, -0.03543202, -0.07060306, 0.10347917)
#'
#' h.prime.temp = NULL #use numerical gradient (optim_grad=FALSE)
#' ForLion_MLM_Optimal(J=J, n.factor=n.factor.temp, factor.level=factor.level.temp, xlist_fix=NULL,
#' hfunc=hfunc.temp,h.prime=h.prime.temp, bvec=bvec.temp, link=link.temp, optim_grad=FALSE)
#'
#'

ForLion_MLM_Optimal <- function(J, n.factor, factor.level,var_names=NULL,xlist_fix=NULL, hfunc, h.prime, bvec, link="continuation", Fi.func=Fi_MLM_func, delta0=1e-5, epsilon=1e-12, reltol=1e-5, delta=0, maxit=100, random=FALSE, nram=3, rowmax=NULL, Xini=NULL, random.initial=FALSE, nram.initial=3, optim_grad=FALSE) {
  d.factor=length(n.factor);             # number of factors
  p.factor=length(bvec);                 # number of predictors
  k.continuous=sum(n.factor==0);         # number of continuous factors
  if(delta==0) delta=reltol;

  #    Case I: all factors are discrete
  if(k.continuous==0) {
    if(is.null(Xini)){ xtemp=xmat_discrete_self(factor.level, rowmax=rowmax)} else {xtemp=Xini;}
    m.design=nrow(xtemp); # initial number of design points
    X.mat = rep(0,J*p.factor*m.design);
    dim(X.mat)=c(J, p.factor, m.design)  # initial model matrix X
    for(i in 1:m.design) {
      if(ncol(xtemp)==1){htemp=hfunc(xtemp[i]);}else{htemp=hfunc(xtemp[i,]);}
      X.mat[,,i]=htemp;
    };

    p0= rep(1/m.design, m.design)
    p0 = p0/sum(p0)
    optemp=liftoneDoptimal_MLM_func(m=m.design, p=p.factor, Xi=X.mat, J=J, thetavec=bvec, link=link, reltol=reltol, maxit=maxit, p00=p0, random=random, nram=nram);

    m.design.ans=m.design=sum(optemp$p>0);       # updated number of design point
    x.factor.ans=x.design=xtemp[optemp$p>0,];  # updated list of design points
    p.ans=p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
    det.ans=det.design=optemp$Maximum;      # optimized determinant
    x.model.ans=X.mat = X.mat[, , optemp$p>0];   # updated model matrix
    convergence.ans=converge.design=optemp$convergence;  # TRUE or FALSE
    itmax.design=0;               # no candidate design points considered
  };                        #     End of Case I

  #    Case II: all factors are continuous
  if(k.continuous==d.factor) {
    lvec=uvec=rep(0, d.factor);     # lower bounds and upper bounds for continuous factors
    for(i in 1:d.factor) {lvec[i]=min(factor.level[[i]]); uvec[i]=max(factor.level[[i]]);};
    if(is.null(Xini)){initial.temp=design_initial_self(k.continuous=k.continuous, factor.level=factor.level, MLM=TRUE,xlist_fix=xlist_fix, lvec=lvec, uvec=uvec, bvec=bvec, link=link, h.func=hfunc, Fi.func=Fi.func, delta0=delta0, epsilon = epsilon, maxit=maxit); xtemp=initial.temp$X; p0=initial.temp$p0} else {xtemp=Xini; p0=NULL}  #no initial design
    if(k.continuous==1) m.design=length(xtemp) else m.design=nrow(xtemp);                   # initial number of design points
    X.mat = rep(0,J*p.factor*m.design);
    dim(X.mat)=c(J, p.factor, m.design)  # initial model matrix X
    for(i in 1:m.design) {
      if(is.null(ncol(xtemp))){htemp=hfunc(xtemp[i]);}else{htemp=hfunc(xtemp[i,]);}
      X.mat[,,i]=htemp;
    };

    if(is.null(p0)){p0 = rep(1/m.design, m.design)}

    optemp=liftoneDoptimal_MLM_func(m=m.design, p=p.factor, Xi=X.mat, J=J, thetavec=bvec, link = link, reltol=reltol, maxit=maxit, p00=p0, random=random, nram=nram)
    # new point step (for each idiscrete, give d and partial d function, use optim to find max d compare to p)
    m.design=sum(optemp$p>0);       # updated number of design point ##deleting w_i = 0
    if(k.continuous==1){x.design=xtemp[optemp$p>0]}else{x.design=xtemp[optemp$p>0,];}  # updated list of design points
    p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
    det.design=optemp$Maximum;      # optimized |X'UX|
    X.mat = X.mat[, ,optemp$p>0];   # updated model matrix

    #calculate F_Xi, Fisher information matrix for current design
    F.mat = 0
    Fi = vector()
    for(i in 1:m.design){
      new_Fi = Fi.func(X.mat[, ,i], bvec, link)
      Fi = append(Fi, list(new_Fi$F_x))
      F.mat = F.mat + p.design[i]*new_Fi$F_x
    }
    inv.F.mat = solve(F.mat, tol=.Machine$double.xmin)

    #calculate d(x, Xi) function
    hfunc1 <- function(y) { hfunc(c(y)); };
    d_x_Xi = function(y){
      hy=hfunc1(y)
      F_newpoint = Fi.func(hy, bvec, link)
      d=psych::tr(inv.F.mat %*% F_newpoint$F_x)
      d=as.numeric(d)
      return(-d)
    }
    grad_d = function(y){
      hy=hfunc1(y)
      Fi_ans = Fi.func(hy, bvec, link)
      Ux = Fi_ans$U_x
      dprime = dprime_func_self(y, bvec, hfunc, h.prime, inv.F.mat, Ux, link=link,k.continuous)
      return(-dprime)
    }
    x0=(lvec+uvec)/2;
    if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
    if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
    #compare with boundary values, if boundary value is better, then use boundary values
    low_dvalue = d_x_Xi(lvec)
    up_dvalue = d_x_Xi(uvec)
    if(low_dvalue < ytemp$value){ytemp$par=lvec; ytemp$value=low_dvalue}
    if(up_dvalue < ytemp$value){ytemp$par=uvec; ytemp$value=up_dvalue}

    #random points
    if(random) for(ia in 1:nram) {
      x0r=x0;
      for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
      if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
      if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
      if(ytemp1$value < ytemp$value) {x0=x0r; ytemp=ytemp1; };
    };
    ystar=c(ytemp$par);
    fvalue=ytemp$value

    nit=1;                        # number of candidate y searched
    while((-ytemp$value/p.factor-1 > reltol)&&(nit < maxit)) {
      # new point considered to be added
      hystar=hfunc(ystar);            # h(ystar)
      alphat=0;                       # calculate weight of ystar

      #add new points to the design update 2023/11/28
      m.design = m.design + 1
      if(k.continuous==1){x.design = c(x.design, ystar)}else{x.design = rbind(x.design, ystar)}
      p.design=c(p.design, alphat);

      #Step 2 merging: calculate distance between design points, if min distance < tolerance merge
      dtemp=as.matrix(stats::dist(x.design));
      diag(dtemp)=Inf;
      atemp=min(dtemp)
      while((atemp<delta)){ # merge closest two neighbors
            #before merging two closest points, save the current state of the design
            x.design_old=x.design
            p.design_old=p.design
            m.design_old=m.design

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
            p.design_mer=c(p.design_old[-c(i1,i2)], p.design_old[i1]+p.design_old[i2])
            m.design_mer=m.design_old-1

            # Build X.mat_mer and calculate the Fisher information matrix (F.mat_mer)
            X.mat_mer = rep(0,J*p.factor*m.design_mer);
            dim(X.mat_mer)=c(J, p.factor, m.design_mer)  # initial model matrix X
            for(i in 1:m.design_mer){
              if(k.continuous==1){htemp_mer=hfunc(x.design_mer[i])}else{htemp_mer=hfunc(x.design_mer[i, ])};
              X.mat_mer[,,i]=htemp_mer;
            };

            #calculate F_Xi, Fisher information matrix for current design
            F.mat_mer = 0
            Fi_mer = vector()
            for(i in 1:m.design_mer){
              new_Fi_mer = Fi.func(X.mat_mer[, ,i], bvec, link)
              Fi_mer = append(Fi_mer, list(new_Fi_mer$F_x))
              F.mat_mer = F.mat_mer + p.design_mer[i]*new_Fi_mer$F_x
            }

            eigen_values<-eigen(F.mat_mer)
            min_engenvalue<-min(eigen_values$values)

            if(min_engenvalue<=epsilon){
              x.design=x.design_old
              p.design=p.design_old
              m.design=m.design_old
              break
            }else{
              x.design=x.design_mer
              p.design=p.design_mer
              m.design=m.design_mer
            }
        dtemp=as.matrix(stats::dist(x.design));
        diag(dtemp)=Inf;
        atemp=min(dtemp)
      }
      X.mat = rep(0,J*p.factor*m.design);
      dim(X.mat)=c(J, p.factor, m.design)  # initial model matrix X
      for(i in 1:m.design){
        if(k.continuous==1){htemp=hfunc(x.design[i])}else{htemp=hfunc(x.design[i, ])};
        X.mat[,,i]=htemp;
      };

      #Go back to step3:lift-one
      optemp=liftoneDoptimal_MLM_func(m=m.design, p=p.factor, Xi=X.mat, J=J, thetavec=bvec,link=link, reltol=reltol, maxit=maxit, p00=p.design, random=random, nram=nram)

      #step4:deleting step new point step
      m.design=sum(optemp$p>0);       # updated number of design point ##deleting w_i = 0
      x.design=if(k.continuous==1){x.design[optemp$p>0]}else{x.design[optemp$p>0,];}  # updated list of design points
      p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
      det.design=optemp$Maximum;      # optimized |X'UX|
      X.mat = X.mat[, ,optemp$p>0];   # updated model matrix

      #step5: new point step (for each idiscrete, give d and partial d function, use optim to find max d compare to p)
      #calculate F_Xi, Fisher information matrix for current design
      F.mat = 0
      Fi = vector()
      for(i in 1:m.design){
        new_Fi = Fi.func(X.mat[, ,i], bvec, link)
        Fi = append(Fi, list(new_Fi$F_x))
        F.mat = F.mat + p.design[i]*new_Fi$F_x
      }
      inv.F.mat = solve(F.mat, tol=.Machine$double.xmin)


      #calculate d(x, Xi) function
      hfunc1 <- function(y) { hfunc(c(y)); };
      d_x_Xi = function(y){
        hy=hfunc1(y)
        F_newpoint = Fi.func(hy, bvec, link)
        d=psych::tr(inv.F.mat %*% F_newpoint$F_x)
        d=as.numeric(d)
        return(-d)
      }
      grad_d = function(y){
        hy=hfunc1(y)
        Fi_ans = Fi.func(hy, bvec, link)
        Ux = Fi_ans$U_x
        dprime = dprime_func_self(y, bvec, hfunc, h.prime, inv.F.mat, Ux, link=link,k.continuous)
        return(-dprime)
      }
      x0=(lvec+uvec)/2;
      if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
      if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
      #compare with boundary values, if boundary value is better, then use boundary values
      low_dvalue = d_x_Xi(lvec)
      up_dvalue = d_x_Xi(uvec)
      if(low_dvalue < ytemp$value){ytemp$par=lvec; ytemp$value=low_dvalue}
      if(up_dvalue < ytemp$value){ytemp$par=uvec; ytemp$value=up_dvalue}

      #random points
      if(random) for(ia in 1:nram) {
        x0r=x0;
        for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
        if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
        if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
        if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
      };

      ystar=c(ytemp$par);
      fvalue=ytemp$value;
      nit=nit+1;                    # number of candidate y searched
    }# end of while(d()>p) loop

    m.design.ans=m.design;      #reported num of design points
    x.factor.ans=x.design #reported design points
    p.ans = p.design #reported design proportions
    det.ans = det.design #reported optimized determinant
    x.model.ans = X.mat #reported model matrix
    itmax.design=nit;
    converge.design=(ytemp$convergence==0);  # TRUE or FALSE
    if(-ytemp$value/p.factor-1 > reltol) converge.design=FALSE;
    convergence.ans=converge.design
    #cat("\n no random initial:\n", "x.factor.ans:", x.factor.ans, "\np.ans:", p.ans, "\ndet.ans:", det.ans)#delete

    #random initial points
    if(random.initial){
      for(num in 1:nram.initial){
        #try different initial points
        initial.temp=design_initial_self(k.continuous=k.continuous, factor.level=factor.level, MLM=TRUE, xlist_fix=xlist_fix, lvec=lvec, uvec=uvec, bvec=bvec, link=link, h.func=hfunc, Fi.func=Fi.func, delta0=delta0, epsilon=epsilon, maxit=maxit); xtemp=initial.temp$X; p0=initial.temp$p0 #random initial design
        if(k.continuous==1) m.design=length(xtemp) else m.design=nrow(xtemp);                   # initial number of design points
        X.mat = rep(0,J*p.factor*m.design);
        dim(X.mat)=c(J, p.factor, m.design)  # initial model matrix X
        for(i in 1:m.design) {
          if(is.null(ncol(xtemp))){htemp=hfunc(xtemp[i]);}else{htemp=hfunc(xtemp[i,]);}
          X.mat[,,i]=htemp;
        };

        optemp=liftoneDoptimal_MLM_func(m=m.design, p=p.factor, Xi=X.mat, J=J, thetavec=bvec, link=link, reltol=reltol, maxit=maxit, p00=p0, random=random, nram=nram)
        # new point step (for each idiscrete, give d and partial d function, use optim to find max d compare to p)
        m.design=sum(optemp$p>0);       # updated number of design point ##deleting w_i = 0
        if(k.continuous==1){x.design=xtemp[optemp$p>0]}else{x.design=xtemp[optemp$p>0,];}  # updated list of design points
        p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
        det.design=optemp$Maximum;      # optimized |X'UX|
        X.mat = X.mat[, ,optemp$p>0];   # updated model matrix

        #calculate F_Xi, Fisher information matrix for current design
        F.mat = 0
        Fi = vector()
        for(i in 1:m.design){
          new_Fi = Fi.func(X.mat[, ,i], bvec, link)
          Fi = append(Fi, list(new_Fi$F_x))
          F.mat = F.mat + p.design[i]*new_Fi$F_x
        }
        inv.F.mat = solve(F.mat, tol=.Machine$double.xmin)

        #calculate d(x, Xi) function
        hfunc1 <- function(y) { hfunc(c(y)); };
        d_x_Xi = function(y){
          hy=hfunc1(y)
          F_newpoint = Fi.func(hy, bvec, link)
          d=psych::tr(inv.F.mat %*% F_newpoint$F_x)
          d=as.numeric(d)
          return(-d)
        }
        grad_d = function(y){
          hy=hfunc1(y)
          Fi_ans = Fi.func(hy, bvec, link)
          Ux = Fi_ans$U_x
          dprime = dprime_func_self(y, bvec, hfunc, h.prime, inv.F.mat, Ux, link=link,k.continuous)
          return(-dprime)
        }
        x0=(lvec+uvec)/2;
        if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
        if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
        #compare with boundary values, if boundary value is better, then use boundary values
        low_dvalue = d_x_Xi(lvec)
        up_dvalue = d_x_Xi(uvec)
        if(low_dvalue < ytemp$value){ytemp$par=lvec; ytemp$value=low_dvalue}
        if(up_dvalue < ytemp$value){ytemp$par=uvec; ytemp$value=up_dvalue}
        #random points
        if(random) for(ia in 1:nram) {
          x0r=x0;
          for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
          if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
          if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
          if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
        };
        ystar=c(ytemp$par);
        fvalue=ytemp$value

        nit=1;                        # number of candidate y searched
        while((-ytemp$value/p.factor-1 > reltol)&&(nit < maxit)) {
          # new point considered to be added
          hystar=hfunc(ystar);            # h(ystar)
          alphat=0;                       # calculate weight of ystar

          #add new points to the design update 2023/11/28
          m.design = m.design + 1
          if(k.continuous==1){x.design = c(x.design, ystar)}else{x.design = rbind(x.design, ystar)}
          p.design=c(p.design, alphat);

          #Step 2 merging: calculate distance between design points, if min distance < tolerance merge
          dtemp=as.matrix(stats::dist(x.design));
          diag(dtemp)=Inf;
          atemp=min(dtemp)
          while((atemp<delta)){ # merge closest two neighbors
            #before merging two closest points, save the current state of the design
            x.design_old=x.design
            p.design_old=p.design
            m.design_old=m.design

            #identify and merge the two closest design points
            i1=which.min(apply(dtemp,1,min)); # index of design point to be merged
            i2=which.min(dtemp[i1,]); # index of design point to be merged
            if(k.continuous==1){
              ystar1=(p.design_old[i1]*x.design_old[i1]+p.design_old[i2]*x.design_old[i2])/(p.design_old[i1]+p.design_old[i2]); # merged design point
              x.design_mer=c(x.design_old[-c(i1,i2)], ystar1)
            }else{
              ystar1=(p.design_old[i1]*x.design_old[i1,]+p.design_old[i2]*x.design_old[i2,])/(p.design_old[i1]+p.design_old[i2]);# merged design point
              x.design_mer=rbind(x.design_old[-c(i1,i2), ], ystar1) # update x.design
            }
            p.design_mer=c(p.design_old[-c(i1,i2)], p.design_old[i1]+p.design_old[i2])
            m.design_mer=m.design_old-1

            # Build X.mat_mer and calculate the Fisher information matrix (F.mat_mer)
            X.mat_mer = rep(0,J*p.factor*m.design_mer);
            dim(X.mat_mer)=c(J, p.factor, m.design_mer)  # initial model matrix X
            for(i in 1:m.design_mer){
              if(k.continuous==1){htemp_mer=hfunc(x.design_mer[i])}else{htemp_mer=hfunc(x.design_mer[i, ])};
              X.mat_mer[,,i]=htemp_mer;
            };

            #calculate F_Xi, Fisher information matrix for current design
            F.mat_mer = 0
            Fi_mer = vector()
            for(i in 1:m.design_mer){
              new_Fi_mer = Fi.func(X.mat_mer[, ,i], bvec, link)
              Fi_mer = append(Fi_mer, list(new_Fi_mer$F_x))
              F.mat_mer = F.mat_mer + p.design_mer[i]*new_Fi_mer$F_x
            }

            eigen_values<-eigen(F.mat_mer)
            min_engenvalue<-min(eigen_values$values)

            if(min_engenvalue<=epsilon){
              x.design=x.design_old
              p.design=p.design_old
              m.design=m.design_old
              break
            }else{
              x.design=x.design_mer
              p.design=p.design_mer
              m.design=m.design_mer
            }
          dtemp=as.matrix(stats::dist(x.design));
          diag(dtemp)=Inf;
          atemp=min(dtemp)
          }
          X.mat = rep(0,J*p.factor*m.design);
          dim(X.mat)=c(J, p.factor, m.design)  # initial model matrix X
          for(i in 1:m.design){
          if(k.continuous==1){htemp=hfunc(x.design[i])}else{htemp=hfunc(x.design[i, ])};
          X.mat[,,i]=htemp;
          };

          #Go back to step3:lift-one
          optemp=liftoneDoptimal_MLM_func(m=m.design,  p=p.factor, Xi=X.mat, J=J, thetavec=bvec,link=link, reltol=reltol, maxit=maxit, p00=p.design, random=random, nram=nram)

          #step4:deleting step new point step
          m.design=sum(optemp$p>0);       # updated number of design point ##deleting w_i = 0
          x.design=if(k.continuous==1){x.design[optemp$p>0]}else{x.design[optemp$p>0,];}  # updated list of design points
          p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
          det.design=optemp$Maximum;      # optimized |X'UX|
          X.mat = X.mat[, ,optemp$p>0];   # updated model matrix
          #cat("\n", num, "th random points:", "\nx.design", x.design, "\np.design",p.design, "\ndet.design", det.design) #delete

          #step5: new point step (for each idiscrete, give d and partial d function, use optim to find max d compare to p)
          #calculate F_Xi, Fisher information matrix for current design
          F.mat = matrix(0, nrow=p.factor, ncol=p.factor)
          Fi = vector()
          for(i in 1:m.design){
            new_Fi = Fi.func(X.mat[, ,i], bvec, link)
            Fi = append(Fi, list(new_Fi$F_x))
            F.mat = F.mat + p.design[i]*new_Fi$F_x
          }
          #cat("\nFmat:", F.mat) #delete
          inv.F.mat = solve(F.mat, tol=.Machine$double.xmin)


          #calculate d(x, Xi) function
          hfunc1 <- function(y) { hfunc(c(y)); };
          d_x_Xi = function(y){
            hy=hfunc1(y)
            F_newpoint = Fi.func(hy, bvec, link)
            d=psych::tr(inv.F.mat %*% F_newpoint$F_x)
            d=as.numeric(d)
            return(-d)
          }
          grad_d = function(y){
            hy=hfunc1(y)
            Fi_ans = Fi.func(hy, bvec, link)
            Ux = Fi_ans$U_x
            dprime = dprime_func_self(y, bvec, hfunc, h.prime, inv.F.mat, Ux, link=link,k.continuous)
            return(-dprime)
          }
          x0=(lvec+uvec)/2;
          if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
          if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
          #compare with boundary values, if boundary value is better, then use boundary values
          low_dvalue = d_x_Xi(lvec)
          up_dvalue = d_x_Xi(uvec)
          if(low_dvalue < ytemp$value){ytemp$par=lvec; ytemp$value=low_dvalue}
          if(up_dvalue < ytemp$value){ytemp$par=uvec; ytemp$value=up_dvalue}
          #random points
          if(random) for(ia in 1:nram) {
            x0r=x0;
            for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
            if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, gr=grad_d,method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
            if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
            if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
          };

          ystar=c(ytemp$par);
          fvalue=ytemp$value;
          nit=nit+1;                    # number of candidate y searched
        }# end of while(d()>p) loop

        #compare the new results with the existing results, replace if better
        if(det.design > det.ans){
          m.design.ans=m.design;      #reported num of design points
          x.factor.ans=x.design #reported design points
          p.ans = p.design #reported design proportions
          det.ans = det.design #reported optimized determinant
          x.model.ans = X.mat #reported model matrix
          itmax.design=nit;
          converge.design=(ytemp$convergence==0);  # TRUE or FALSE
          if(-ytemp$value/p.factor-1 > reltol) converge.design=FALSE;
          convergence.ans=converge.design
        } # end of if condition (change better random points)

      }#end of for loop nram.initial

    } #end of if(random.initial)
  }# end of Case II

  #    Case III: some factors are continuous #firstly change Case III, then generalize to Case I and Case II
  if((k.continuous>0)&&(k.continuous<d.factor)) {
    lvec=uvec=rep(0, k.continuous);     # lower bounds and upper bounds for continuous factors
    for(i in 1:k.continuous) {lvec[i]=min(factor.level[[i]]); uvec[i]=max(factor.level[[i]]);}; #read in continuous covariates boundary
    if(is.null(Xini)){initial.temp=design_initial_self(k.continuous=k.continuous, factor.level=factor.level, MLM=TRUE, xlist_fix=xlist_fix, lvec=lvec, uvec=uvec, bvec=bvec, link=link, h.func=hfunc, Fi.func=Fi.func, delta0=delta0, epsilon = epsilon, maxit=500); xtemp=initial.temp$X; p0=initial.temp$p0} else {xtemp=Xini; p0=NULL}  #no initial design

    ## update on 2022/08/28 change Xw.discrete.self function to sequentially and randomly choose x_i^0 => design_initial_self function

    m.design=nrow(xtemp) # initial number of design points
    X.mat = rep(0,J*p.factor*m.design);
    dim(X.mat)=c(J, p.factor, m.design)  # initial model matrix X
    for(i in 1:m.design) {
      if(ncol(xtemp)==1){htemp=hfunc(xtemp[i]);}else{htemp=hfunc(xtemp[i,]);}
      X.mat[,,i]=htemp;
    };
    if(is.null(p0)){p0 = rep(1/m.design, m.design)}

    optemp=liftoneDoptimal_MLM_func(m=m.design, p=p.factor, Xi=X.mat, J=J, thetavec=bvec, link=link, reltol=reltol, maxit=maxit, p00=p0, random=random, nram=nram)

    # new point step (for each idiscrete, give d and partial d function, use optim to find max d compare to p)
    m.design=sum(optemp$p>0);       # updated number of design point ##deleting w_i = 0
    x.design=xtemp[optemp$p>0,]; # updated list of design points
    p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
    det.design=optemp$Maximum;      # optimized |X'UX|
    X.mat = X.mat[, ,optemp$p>0];   # updated model matrix

    #calculate F_Xi, Fisher information matrix for current design
    F.mat = 0
    Fi = vector()
    for(i in 1:m.design){
      new_Fi = Fi.func(X.mat[, ,i], bvec, link)
      Fi = append(Fi, list(new_Fi$F_x))
      F.mat = F.mat + p.design[i]*new_Fi$F_x
    }
    inv.F.mat = solve(F.mat, tol=.Machine$double.xmin)

    #calculate d(x, Xi) function
    if(is.null(xlist_fix)){xdiscrete=xmat_discrete_self(factor.level[(k.continuous+1):d.factor]);}else{xdiscrete=xlist_fix;}
    ndiscrete=dim(xdiscrete)[1];

    for(idiscrete in 1:ndiscrete) {
      hfunc1 <- function(y) { hfunc(c(y, xdiscrete[idiscrete,])); };
      d_x_Xi = function(y){
        hy=hfunc1(y)
        F_newpoint = Fi.func(hy, bvec, link)
        d=psych::tr(inv.F.mat %*% F_newpoint$F_x)
        d=as.numeric(d)
        return(-d)
      }
      grad_d = function(y){
        hy=hfunc1(y)
        Fi_ans = Fi.func(hy, bvec, link)
        Ux = Fi_ans$U_x
        dprime = dprime_func_self(c(y, xdiscrete[idiscrete,]), bvec, hfunc, h.prime, inv.F.mat, Ux, link=link,k.continuous)
        return(-dprime)
      }
      x0=(lvec+uvec)/2;
      if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
      if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
      # ytemp = spg(par=x0, fn=d_x_Xi, lower=lvec, upper=uvec, control=list(maxit=maxit))
      #compare with boundary values, if boundary value is better, then use boundary values
      low_dvalue = d_x_Xi(lvec)
      up_dvalue = d_x_Xi(uvec)
      if(low_dvalue < ytemp$value){ytemp$par=lvec; ytemp$value=low_dvalue}
      if(up_dvalue < ytemp$value){ytemp$par=uvec; ytemp$value=up_dvalue}
      ytempstar=ytemp;
      #random points
      if(random) for(ia in 1:nram) {
        x0r=x0;
        for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
        if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
        if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
        # ytemp1=spg(par=x0r, fn=d_x_Xi, lower=lvec, upper=uvec, control=list(maxit=maxit));
        if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
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

    } #end of for loop of idiscrete


    nit = 1
    while((-fvalue/p.factor-1 > reltol)&&(nit < maxit)) {
      # new point considered to be added
      hystar=hfunc(ystar);            # h(ystar)
      alphat=0;                       # calculate weight of ystar

      #add new points to the design update 2023/11/28
      m.design = m.design + 1
      x.design = rbind(x.design, ystar)
      p.design=c(p.design, alphat);

      # Step 2 merging: calculate distance between design points, if min distance < tolerance
      dtemp=as.matrix(stats::dist(x.design));
      diag(dtemp)=Inf;
      atemp=min(dtemp)
      while((atemp<delta)){ # merge closest two neighbors
            #before merging two closest points, save the current state of the design
            x.design_old=x.design
            p.design_old=p.design
            m.design_old=m.design

            #identify and merge the two closest design points
            i1=which.min(apply(dtemp,1,min)); # index of design point to be merged
            i2=which.min(dtemp[i1,]); # index of design point to be merged
            ystar1=(p.design_old[i1]*x.design_old[i1,]+p.design_old[i2]*x.design_old[i2,])/(p.design_old[i1]+p.design_old[i2]); # merged design point
            x.design_mer=rbind(x.design_old[-c(i1,i2),], ystar1); # update x.design
            p.design_mer=c(p.design_old[-c(i1,i2)], p.design_old[i1]+p.design_old[i2])
            m.design_mer=m.design_old-1

            # Build X.mat_mer and calculate the Fisher information matrix (F.mat_mer)
            X.mat_mer = rep(0,J*p.factor*m.design_mer);
            dim(X.mat_mer)=c(J, p.factor, m.design_mer)  # initial model matrix X
            for(i in 1:m.design_mer){
              htemp_mer=hfunc(x.design_mer[i, ]);
              X.mat_mer[,,i]=htemp_mer;
            };

            #calculate F_Xi, Fisher information matrix for current design
            F.mat_mer = 0
            Fi_mer = vector()
            for(i in 1:m.design_mer){
              new_Fi_mer = Fi.func(X.mat_mer[, ,i], bvec, link)
              Fi_mer = append(Fi_mer, list(new_Fi_mer$F_x))
              F.mat_mer = F.mat_mer + p.design_mer[i]*new_Fi_mer$F_x
            }

            eigen_values<-eigen(F.mat_mer)
            min_engenvalue<-min(eigen_values$values)

            if(min_engenvalue<=epsilon){
              x.design=x.design_old
              p.design=p.design_old
              m.design=m.design_old
              break
            }else{
              x.design=x.design_mer
              p.design=p.design_mer
              m.design=m.design_mer
            }
        dtemp=as.matrix(stats::dist(x.design));
        diag(dtemp)=Inf;
        atemp=min(dtemp)
      }
      X.mat = rep(0,J*p.factor*m.design);
      dim(X.mat)=c(J, p.factor, m.design)  # initial model matrix X
      for(i in 1:m.design){
        htemp=hfunc(x.design[i,]);
        X.mat[,,i]=htemp;
      };

      #Go back to step3:lift-one
      optemp=liftoneDoptimal_MLM_func(m=m.design, p=p.factor, Xi=X.mat, J=J, thetavec=bvec, link=link, reltol=reltol, maxit=maxit, p00=p.design, random=random, nram=nram)

      #step4:deleting step new point step
      m.design=sum(optemp$p>0);       # updated number of design point ##deleting w_i = 0
      x.design=x.design[optemp$p>0,];  # updated list of design points
      p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
      det.design=optemp$Maximum;      # optimized |X'UX|
      X.mat = X.mat[, ,optemp$p>0];   # updated model matrix

      #step5: new point step (for each idiscrete, give d and partial d function, use optim to find max d compare to p)
      #calculate F_Xi, Fisher information matrix for current design
      F.mat = 0
      Fi = vector()
      for(i in 1:m.design){
        new_Fi = Fi.func(X.mat[, ,i], bvec, link)
        Fi = append(Fi, list(new_Fi$F_x))
        F.mat = F.mat + p.design[i]*new_Fi$F_x
      }
      inv.F.mat = solve(F.mat, tol=.Machine$double.xmin)

      #calculate d(x, Xi) function
      if(is.null(xlist_fix)){xdiscrete=xmat_discrete_self(factor.level[(k.continuous+1):d.factor]);}else{xdiscrete=xlist_fix;}
      #xdiscrete=xmat_discrete_self(factor.level[(k.continuous+1):d.factor]);
      ndiscrete=dim(xdiscrete)[1];
      for(idiscrete in 1:ndiscrete) {
        hfunc1 <- function(y) { hfunc(c(y, xdiscrete[idiscrete,])); };
        d_x_Xi = function(y){
          hy=hfunc1(y)
          F_newpoint = Fi.func(hy, bvec, link)
          d=psych::tr(inv.F.mat %*% F_newpoint$F_x)
          d=as.numeric(d)
          return(-d)
        }
        grad_d = function(y){
          hy=hfunc1(y)
          Fi_ans = Fi.func(hy, bvec, link)
          Ux = Fi_ans$U_x
          dprime = dprime_func_self(c(y, xdiscrete[idiscrete,]), bvec, hfunc, h.prime, inv.F.mat, Ux, link=link,k.continuous)
          return(-dprime)
        }
        x0=(lvec+uvec)/2;
        if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
        if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
        # ytemp=spg(par=x0, fn=d_x_Xi, lower=lvec, upper=uvec, control=list(maxit=maxit));

        #compare with boundary values, if boundary value is better, then use boundary values
        low_dvalue = d_x_Xi(lvec)
        up_dvalue = d_x_Xi(uvec)
        if(low_dvalue < ytemp$value){ytemp$par=lvec; ytemp$value=low_dvalue}
        if(up_dvalue < ytemp$value){ytemp$par=uvec; ytemp$value=up_dvalue}
        #random points
        if(random) for(ia in 1:nram) {
          x0r=x0;
          for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
          if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
          if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
          # ytemp1=spg(par=x0r, fn=d_x_Xi, lower=lvec, upper=uvec, control=list(maxit=maxit));
          if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
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

      } #end of for loop of idiscrete

      nit=nit+1;                    # number of candidate y searched
    }# end of while(d()>p) loop

    itmax.design=nit;
    converge.design=(ytempstar$convergence==0);  # TRUE or FALSE
    if(-ytempstar$value/p.factor-1 > reltol) converge.design=FALSE;
    #updated on 10/23/2022 record as ....ans and to compare with results from random initial points
    m.design.ans=m.design;      #reported num of design points
    x.factor.ans=x.design #reported design points
    p.ans = p.design #reported design proportions
    det.ans = det.design #reported optimized determinant
    x.model.ans = X.mat #reported model matrix
    itmax.design=nit;
    converge.design=(ytempstar$convergence==0);  # TRUE or FALSE
    if(-ytempstar$value/p.factor-1 > reltol) converge.design=FALSE;
    convergence.ans = converge.design

    #update on 08/30/2022 add random initial x_i^(0)
    if(random.initial){
      for(num in 1:nram.initial){
        #try different random x_i^0
        initial.temp=design_initial_self(k.continuous=k.continuous, factor.level=factor.level, MLM=TRUE, xlist_fix=xlist_fix, lvec=lvec, uvec=uvec, bvec=bvec, link=link, h.func=hfunc, Fi.func=Fi.func, delta0=delta0, epsilon = epsilon, maxit=500); xtemp=initial.temp$X; p0=initial.temp$p0  #random initial design

        ## update on 2022/08/28 change Xw.discrete.self function to sequentially and randomly choose x_i^0 => design_initial_self function

        m.design=ifelse(is.null(nrow(xtemp)), length(xtemp), nrow(xtemp)); # initial number of design points
        X.mat = rep(0,J*p.factor*m.design);
        dim(X.mat)=c(J, p.factor, m.design)  # initial model matrix X
        for(i in 1:m.design) {
          if(ncol(xtemp)==1){htemp=hfunc(xtemp[i]);}else{htemp=hfunc(xtemp[i,]);}
          X.mat[,,i]=htemp;
        };

        optemp=liftoneDoptimal_MLM_func(m=m.design, p=p.factor, Xi=X.mat, J=J, thetavec=bvec, link=link, reltol=reltol, maxit=maxit, p00=p0, random=random, nram=nram)

        # new point step (for each idiscrete, give d and partial d function, use optim to find max d compare to p)
        m.design=sum(optemp$p>0);       # updated number of design point ##deleting w_i = 0
        x.design=xtemp[optemp$p>0,];  # updated list of design points
        p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
        det.design=optemp$Maximum;      # optimized |X'UX|
        X.mat = X.mat[, ,optemp$p>0];   # updated model matrix

        #calculate F_Xi, Fisher information matrix for current design
        F.mat = 0
        Fi = vector()
        for(i in 1:m.design){
          new_Fi = Fi.func(X.mat[, ,i], bvec, link)
          Fi = append(Fi, list(new_Fi$F_x))
          F.mat = F.mat + p.design[i]*new_Fi$F_x
        }
        inv.F.mat = solve(F.mat, tol=.Machine$double.xmin)

        #calculate d(x, Xi) function
        if(is.null(xlist_fix)){xdiscrete=xmat_discrete_self(factor.level[(k.continuous+1):d.factor]);}else{xdiscrete=xlist_fix;}
        #xdiscrete=xmat_discrete_self(factor.level[(k.continuous+1):d.factor]);
        ndiscrete=dim(xdiscrete)[1];
        for(idiscrete in 1:ndiscrete) {
          hfunc1 <- function(y) { hfunc(c(y, xdiscrete[idiscrete,])); };
          d_x_Xi = function(y){
            hy=hfunc1(y)
            F_newpoint = Fi.func(hy, bvec, link)
            d=psych::tr(inv.F.mat %*% F_newpoint$F_x)
            d=as.numeric(d)
            return(-d)
          }
          grad_d = function(y){
            hy=hfunc1(y)
            Fi_ans = Fi.func(hy, bvec, link)
            Ux = Fi_ans$U_x
            dprime = dprime_func_self(c(y, xdiscrete[idiscrete,]), bvec, hfunc, h.prime, inv.F.mat, Ux, link=link,k.continuous)
            return(-dprime)
          }
          x0=(lvec+uvec)/2;
          if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
          if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
          # ytemp=spg(par=x0, fn=d_x_Xi, lower=lvec, upper=uvec, control=list(maxit=maxit));
          #compare with boundary values, if boundary value is better, then use boundary values
          low_dvalue = d_x_Xi(lvec)
          up_dvalue = d_x_Xi(uvec)
          if(low_dvalue < ytemp$value){ytemp$par=lvec; ytemp$value=low_dvalue}
          if(up_dvalue < ytemp$value){ytemp$par=uvec; ytemp$value=up_dvalue}
          #random points
          if(random) for(ia in 1:nram) {
            x0r=x0;
            for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
            if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
            if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
            # ytemp1=spg(par=x0r, fn=d_x_Xi, lower=lvec, upper=uvec, control=list(maxit=maxit));
            if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
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

        } #end of for loop of idiscrete

        nit = 1
        while((-fvalue/p.factor-1 > reltol)&&(nit < maxit)) {
          # add new point into the design
          hystar=hfunc(ystar);            # h(ystar)
          alphat=0;                       # calculate weight of ystar

          #add new points to the design update 2023/11/28
          m.design = m.design + 1
          x.design = rbind(x.design, ystar)
          p.design=c(p.design, alphat);

          #Step 2 merging: calculate distance between design points, if min distance < tolerance
          dtemp=as.matrix(stats::dist(x.design));
          diag(dtemp)=Inf;
          atemp=min(dtemp)
          while((atemp<delta)){ # merge closest two neighbors
            #before merging two closest points, save the current state of the design
            x.design_old=x.design
            p.design_old=p.design
            m.design_old=m.design

            #identify and merge the two closest design points
            i1=which.min(apply(dtemp,1,min)); # index of design point to be merged
            i2=which.min(dtemp[i1,]); # index of design point to be merged
            ystar1=(p.design_old[i1]*x.design_old[i1,]+p.design_old[i2]*x.design_old[i2,])/(p.design_old[i1]+p.design_old[i2]); # merged design point
            x.design_mer=rbind(x.design_old[-c(i1,i2),], ystar1); # update x.design
            p.design_mer=c(p.design_old[-c(i1,i2)], p.design_old[i1]+p.design_old[i2])
            m.design_mer=m.design_old-1

            # Build X.mat_mer and calculate the Fisher information matrix (F.mat_mer)
            X.mat_mer = rep(0,J*p.factor*m.design_mer);
            dim(X.mat_mer)=c(J, p.factor, m.design_mer)  # initial model matrix X
            for(i in 1:m.design_mer){
              htemp_mer=hfunc(x.design_mer[i, ]);
              X.mat_mer[,,i]=htemp_mer;
            };

            #calculate F_Xi, Fisher information matrix for current design
            F.mat_mer = 0
            Fi_mer = vector()
            for(i in 1:m.design_mer){
              new_Fi_mer = Fi.func(X.mat_mer[, ,i], bvec, link)
              Fi_mer = append(Fi_mer, list(new_Fi_mer$F_x))
              F.mat_mer = F.mat_mer + p.design_mer[i]*new_Fi_mer$F_x
            }

            eigen_values<-eigen(F.mat_mer)
            min_engenvalue<-min(eigen_values$values)

            if(min_engenvalue<=epsilon){
              x.design=x.design_old
              p.design=p.design_old
              m.design=m.design_old
              break
            }else{
              x.design=x.design_mer
              p.design=p.design_mer
              m.design=m.design_mer
            }
        dtemp=as.matrix(stats::dist(x.design));
        diag(dtemp)=Inf;
        atemp=min(dtemp)
      }
          X.mat = rep(0,J*p.factor*m.design);
          dim(X.mat)=c(J, p.factor, m.design)  # initial model matrix X
          for(i in 1:m.design){
            htemp=hfunc(x.design[i,]);
            X.mat[,,i]=htemp;
          };

          #Go back to step3:lift-one
          optemp=liftoneDoptimal_MLM_func(m=m.design, p=p.factor, Xi=X.mat, J=J, thetavec=bvec, link=link, reltol=reltol, maxit=maxit, p00=p.design, random=random, nram=nram)

          #step4:deleting step new point step
          m.design=sum(optemp$p>0);       # updated number of design point ##deleting w_i = 0
          x.design=x.design[optemp$p>0,];  # updated list of design points
          p.design=optemp$p[optemp$p>0];  # optimal allocation on current design points
          det.design=optemp$Maximum;      # optimized |X'UX|
          X.mat = X.mat[, ,optemp$p>0];   # updated model matrix

          #step5: new point step (for each idiscrete, give d and partial d function, use optim to find max d compare to p)
          #calculate F_Xi, Fisher information matrix for current design
          F.mat = 0
          Fi = vector()
          for(i in 1:m.design){
            new_Fi = Fi.func(X.mat[, ,i], bvec, link)
            Fi = append(Fi, list(new_Fi$F_x))
            F.mat = F.mat + p.design[i]*new_Fi$F_x
          }
          inv.F.mat = solve(F.mat, tol=.Machine$double.xmin)

          #calculate d(x, Xi) function
          if(is.null(xlist_fix)){xdiscrete=xmat_discrete_self(factor.level[(k.continuous+1):d.factor]);}else{xdiscrete=xlist_fix;}
          # xdiscrete=xmat_discrete_self(factor.level[(k.continuous+1):d.factor]);
          ndiscrete=dim(xdiscrete)[1];
          for(idiscrete in 1:ndiscrete) {
            hfunc1 <- function(y) { hfunc(c(y, xdiscrete[idiscrete,])); };
            d_x_Xi = function(y){
              hy=hfunc1(y)
              F_newpoint = Fi.func(hy, bvec, link)
              d=psych::tr(inv.F.mat %*% F_newpoint$F_x)
              d=as.numeric(d)
              return(-d)
            }
            grad_d = function(y){
              hy=hfunc1(y)
              Fi_ans = Fi.func(hy, bvec, link)
              Ux = Fi_ans$U_x
              dprime = dprime_func_self(c(y, xdiscrete[idiscrete,]), bvec, hfunc, h.prime, inv.F.mat, Ux, link=link,k.continuous)
              return(-dprime)
            }
            x0=(lvec+uvec)/2;
            if(optim_grad==TRUE){ytemp=stats::optim(par=x0, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
            if(optim_grad==FALSE){ytemp=stats::optim(par=x0, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
            # ytemp=spg(par=x0, fn=d_x_Xi, lower=lvec, upper=uvec, control=list(maxit=maxit));
            #compare with boundary values, if boundary value is better, then use boundary values
            low_dvalue = d_x_Xi(lvec)
            up_dvalue = d_x_Xi(uvec)
            if(low_dvalue < ytemp$value){ytemp$par=lvec; ytemp$value=low_dvalue}
            if(up_dvalue < ytemp$value){ytemp$par=uvec; ytemp$value=up_dvalue}
            #random points
            if(random) for(ia in 1:nram) {
              x0r=x0;
              for(i in 1:k.continuous) x0r[i]=lvec[i]+stats::rbeta(1, 0.5, 0.5)*(uvec[i]-lvec[i]);
              if(optim_grad==TRUE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, gr=grad_d, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
              if(optim_grad==FALSE){ytemp1=stats::optim(par=x0r, fn=d_x_Xi, method="L-BFGS-B", lower=lvec, upper=uvec, control=list(maxit=maxit, factr=reltol*1e13));}
              # ytemp1=spg(par=x0r, fn=d_x_Xi, lower=lvec, upper=uvec, control=list(maxit=maxit));
              if(ytemp1$value < ytemp$value) { x0=x0r; ytemp=ytemp1; };
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

          } #end of for loop of idiscrete

          nit=nit+1;                    # number of candidate y searched
        }# end of while(d()>p) loop

        itmax.design=nit;
        converge.design=(ytempstar$convergence==0);  # TRUE or FALSE
        if(-ytempstar$value/p.factor-1 > reltol) converge.design=FALSE;
        convergence.ans=converge.design
        #compare the new results with the existing results, replace if better
        if(det.design > det.ans){
          m.design.ans=m.design;      #reported num of design points
          x.factor.ans=x.design #reported design points
          p.ans = p.design #reported design proportions
          det.ans = det.design #reported optimized determinant
          x.model.ans = X.mat #reported model matrix
          itmax.design=nit;
          converge.design=(ytempstar$convergence==0);  # TRUE or FALSE
          if(-ytempstar$value/p.factor-1 > reltol) converge.design=FALSE;
          convergence.ans=converge.design
        } # end of if condition (change better random points)

      }# end of for loop of nram.initial

    }# end of if(random.initial loop)

  };                                # end of Case III

  #final summarization #updated on 08/30/2022 change name to name.ans
  rownames(x.factor.ans)=NULL;
  rownames(x.model.ans)=NULL;
  dtemp=as.matrix(stats::dist(x.factor.ans)); #merge step
  diag(dtemp)=Inf;
  min.diff=min(dtemp);
  i1=which.min(apply(dtemp,1,min));
  i2=which.min(dtemp[i1,]);
  if(d.factor==1) x.close=x.factor.ans[c(i1,i2)] else {
    x.close=x.factor.ans[c(i1,i2),];
  };
  # x.factor.ans = x.factor.ans[order(x.factor.ans[,1],decreasing=FALSE),]
  # p.ans = p.ans[order(x.factor.ans[,1],decreasing=FALSE)]
  # x.model.ans = x.model.ans[,,order(x.factor.ans[,1],decreasing=FALSE)]
  #cat("\nconverge.ans:", convergence.ans, "\ncalculation:", -ytemp$value/p.factor-1)#delete
  # list(m=m.design.ans, x.factor=x.factor.ans, p=p.ans, det=det.ans, x.model=x.model.ans, convergence=convergence.ans, min.diff=min.diff, x.close=x.close, itmax.design=itmax.design); #updated on 08/30/2022 change name to name.ans

  #define S3 class
  output<-list(m=m.design.ans, x.factor=x.factor.ans, p=p.ans, var.names=var_names, det=det.ans, convergence=convergence.ans, min.diff=min.diff, x.close=x.close, itmax=itmax.design); #updated on 08/30/2022 change name to name.ans
  class(output) <- "design_output"
  return(output)
}
