#' function to generate a initial EW Design for generalized linear models
#'
#' @param k.continuous number of continuous variables
#' @param factor.level lower, upper limit of continuous variables, and discrete levels of categorical variables, continuous factors come first
#' @param Integral_based TRUE or FALSE, if TRUE then we will find the integral-based EW D-optimality otherwise we will find the sample-based EW D-optimality
#' @param b_matrix     The matrix of the sampled parameter values of beta
#' @param joint_Func_b The prior joint probability distribution of the parameters
#' @param Lowerbounds The lower limit of the prior distribution for each parameter
#' @param Upperbounds The upper limit of the prior distribution for each parameter
#' @param xlist_fix the restricted discrete settings to be chosen, default to NULL, if NULL, will generate a discrete uniform random variables
#' @param lvec lower limit of continuous variables
#' @param uvec upper limit of continuous variables
#' @param h.func function, is used to transfer the design point to model matrix (e.g. add interaction term, add intercept)
#' @param link link function, default "continuation", other options "baseline", "adjacent" and "cumulative"
#' @param delta tuning parameter, the distance threshold, || x_i(0) - x_j(0) || >= delta
#' @param epsilon determining f.det > 0 numerically, f.det <= epsilon will be considered as f.det <= 0
#' @param maxit maximum number of iterations
#'
#' @return X        matrix of initial design point
#' @return p0       initial random approximate allocation
#' @return f.det    the determinant of the expected Fisher information matrix for the initial design
#' @export
#'

EW_design_initial_GLM<- function(k.continuous, factor.level, Integral_based, b_matrix, joint_Func_b,Lowerbounds, Upperbounds,xlist_fix=NULL, lvec, uvec, h.func,link="continuation", delta=1e-6, epsilon=1e-12, maxit=1000){
  d.rv = length(factor.level) #number of variables
  if(k.continuous > 0 && (d.rv-k.continuous) > 0){ #mixed case
    #generate initial continuous uniform r.v for continuous variables
    continuous.var = stats::runif(k.continuous, min=lvec, max=uvec)
    #generate initial discrete uniform r.v. for categorical variables
    if(is.null(xlist_fix)){categorical.var = discrete_rv_self(d.rv-k.continuous, factor.level[k.continuous+1:d.rv])}
    else{categorical.var = xlist_fix[sample(nrow(xlist_fix),size=1,replace=TRUE),]}
    x = c(continuous.var, categorical.var) #combine the initial point so it has d variables
  }
  # if(k.continuous == 0){ #discrete case
  #   #generate initial discrete uniform r.v. for categorical variables
  #   categorical.var = discrete_rv_self(d.rv-k.continuous, factor.level[k.continuous+1:d.rv])
  #   x = categorical.var
  # }
  if(d.rv==k.continuous){ #continuous case
    #generate initial continuous uniform r.v for continuous variables
    continuous.var = stats::runif(k.continuous, min=lvec, max=uvec)
    x = continuous.var
  }

  #calculate x's fisher determinant at the beginning (first point)
  m0 = 1 #number of design points currently in the model matrix
  f.det = 0
  iter = 0
  #while loop, check whether det(F) > 0, if not generate new point, record m0
  while(!(f.det > epsilon && m0>2) && (iter<= maxit)){
    if(k.continuous > 0 && (d.rv-k.continuous) > 0){ #mixed case
      #generate initial continuous uniform r.v for continuous variables
      continuous.var = stats::runif(k.continuous, min=lvec, max=uvec)
      #generate initial discrete uniform r.v. for categorical variables
      if(is.null(xlist_fix)){categorical.var = discrete_rv_self(d.rv-k.continuous, factor.level[k.continuous+1:d.rv])}
      else{categorical.var = xlist_fix[sample(nrow(xlist_fix),size=1,replace=TRUE),]}
      new.point = c(continuous.var, categorical.var) #combine the initial point so it has d variables
    }
    # if(k.continuous == 0){ #discrete case
    #   #generate initial discrete uniform r.v. for categorical variables
    #   categorical.var = discrete_rv_self(d.rv-k.continuous, factor.level[k.continuous+1:d.rv])
    #   new.point = categorical.var
    # }
    if(d.rv==k.continuous){ #continuous case
      #generate initial continuous uniform r.v for continuous variables
      continuous.var = stats::runif(k.continuous, min=lvec, max=uvec)
      new.point = continuous.var
    }
    dist = rep(0, m0)
    # for(j in 1:m0){
    #   dist[j] = sqrt(sum((x[j] - new.point)^2))
    # }
    for(j in 1:m0){
      if(d.rv==1){dist[j] = sqrt(sum((x[j] - new.point)^2))}
      else{
        if(m0==1){dist[j]=sqrt(sum((x - new.point)^2))}
        else{
          dist[j] = sqrt(sum((x[j, ] - new.point)^2))}
      }
    }

    #if new point meets the requirements, append to x matrix, update f.det; if not skip
    if(sum(dist < delta)==0){
      m0 = m0+1
      if(d.rv==1){x = c(x, new.point)}else{x = rbind(x, new.point)}
      #x = rbind(x, new.point) #add new row of design points
      p0 = rep(1/m0,m0) #generate approx design as exp r.v.
      #update f.det
      if(Integral_based==TRUE){
        if(d.rv==1){m.design=length(x)} else {m.design=nrow(x);}# initial number of design points
        p.factor=length(Lowerbounds)
        X.mat = matrix(0, m.design, p.factor);  # initial model matrix X
        E_w.vec = rep(0, m.design);     # E_w vector
        for(i in 1:m.design) {
          if(d.rv==1) htemp=EW_Xw_maineffects_self(x=x[i],Integral_based=Integral_based,joint_Func_b=joint_Func_b,Lowerbounds=Lowerbounds, Upperbounds=Upperbounds,link=link, h.func=h.func) else {
            htemp=EW_Xw_maineffects_self(x=x[i,],Integral_based=Integral_based,joint_Func_b=joint_Func_b, Lowerbounds=Lowerbounds, Upperbounds=Upperbounds, link=link, h.func=h.func);
          };
          X.mat[i,]=htemp$X;
          E_w.vec[i]=htemp$E_w;
        };
        f.det = det(t(X.mat * (p0*E_w.vec)) %*% X.mat)
      }else{
        if(d.rv==1){m.design=length(x)} else {m.design=nrow(x);}# initial number of design points
        p.factor=dim(b_matrix)[2]
        X.mat = matrix(0, m.design, p.factor);  # initial model matrix X
        E_w.vec = rep(0, m.design);     # E_w vector
        for(i in 1:m.design) {
          if(d.rv==1)   htemp=EW_Xw_maineffects_self(x=x[i],Integral_based=Integral_based,b_matrix=b_matrix,link=link, h.func=h.func) else {
            htemp=EW_Xw_maineffects_self(x=x[i,],Integral_based=Integral_based,b_matrix=b_matrix, link=link, h.func=h.func);
          };
          X.mat[i,]=htemp$X;
          E_w.vec[i]=htemp$E_w;
        };
        f.det=det(t(X.mat * (p0*E_w.vec)) %*% X.mat)    # X^T W X
      }
    } #end of if
    iter = iter + 1
  } #end of while loop

  if(d.rv==1){X = x[do.call(order, as.data.frame(x))]}else
  {X = x[do.call(order, as.data.frame(x)), ]}
  #generate approx design as exp r.v.
  #exp.rv = rexp(m0)
  #p0 = exp.rv / sum(exp.rv)

  list(X = X, p0 = p0, f.det=f.det)
}#end of function
