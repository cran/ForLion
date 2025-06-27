#' function to generate initial design with design points and the approximate allocation
#' @param k.continuous number of continuous variables
#' @param factor.level list of distinct factor levels, “(min, max)” for continuous factors that always come first, finite sets for discrete factors.
#' @param MLM TRUE or FALSE, TRUE: generate initial design for multinomial logistic model, FALSE: generate initial design for generalized linear model
#' @param xlist_fix list of discrete factor experimental settings under consideration, default NULL indicating a list of all possible discrete factor experimental settings will be used.
#' @param lvec lower limit of continuous variables
#' @param uvec upper limit of continuous variables
#' @param bvec assumed parameter values of beta
#' @param link link function, default "continuation", other options "baseline", "adjacent" and "cumulative"
#' @param h.func function for generating the corresponding model matrix or predictor vector, given an experimental setting or design point.
#' @param Fi.func function, is used to calculate Fisher inforamtion for a design point, default to be Fi_MLM_func() in the package
#' @param delta0 tuning parameter, the distance threshold, || x_i(0) - x_j(0) || >= delta0
#' @param epsilon tuning parameter as converging threshold, such that, a nonnegative number is regarded as numerical zero if less than epsilon, default 1e-12.
#' @param maxit maximum number of iterations
#'
#' @return X      matrix of initial design point
#' @return p0     initial random approximate allocation
#' @return f.det  the determinant of Fisher information matrix for the initial design
#' @export
#'
#' @examples
#' k.continuous.temp=5
#' link.temp = "cumulative"
#' n.factor.temp = c(0,0,0,0,0,2)  # 1 discrete factor w/ 2 levels + 5 continuous
#' ## Note: Always put continuous factors ahead of discrete factors,
#' ## pay attention to the order of coefficients paring with predictors
#' lvec.temp = c(-25,-200,-150,-100,0,-1)
#' uvec.temp = c(25,200,0,0,16,1)
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
#' design_initial_self(k.continuous=k.continuous.temp, factor.level=n.factor.temp,
#' MLM=TRUE,xlist_fix=NULL, lvec=lvec.temp,uvec=uvec.temp, bvec=bvec.temp,
#' h.func=hfunc.temp,link=link.temp)
#'
#'


design_initial_self = function(k.continuous, factor.level, MLM, xlist_fix=NULL, lvec, uvec, bvec, h.func, link="continuation", Fi.func=Fi_MLM_func, delta0=1e-6, epsilon=1e-12, maxit=1000){
  ## function to generate initial design x_i^(0) and w_i^(0) #update 2022/08/28
  ## input:
  ##        k.continuous: number of continuous variables
  ##        factor.level: lower, upper limit of continuous variables, and discrete levels of categorical variables
  ##        lvec: lower limit of continuous variables
  ##        uvec: upper limit of continuous variables
  ##        bvec: assumed parameter values of beta, same length of h(y)
  ##        link: link function, default "continuation"
  ##        h.func: function, is used to transfer the design matrix to model matrix (e.g. add interaction term, add intercept)
  ##        Fi.func: function, is used to calculate rowwise fisher information Fi
  ##        delta0: tuning parameter, the distance threshold, || x_i(0) - x_j(0) || >= delta0
  ##        epsilon: for determining f.det > 0 numerically, f.det <= epsilon will be considered as f.det <= 0
  ## output:
  ##        X: matrix of initial design points
  ##        p0: initial approx design of X
  d.rv = length(factor.level) #number of variables
  if(k.continuous > 0 && (d.rv-k.continuous) > 0){ #mixed case
    #generate initial continuous uniform r.v for continuous variables
    continuous.var = stats::runif(k.continuous, min=lvec, max=uvec)
    #generate initial discrete uniform r.v. for categorical variables
    if(is.null(xlist_fix)){categorical.var = discrete_rv_self(d.rv-k.continuous, factor.level[k.continuous+1:d.rv])}
    else{categorical.var = xlist_fix[sample(nrow(xlist_fix),size=1,replace=TRUE),]} ##only apply for more than 1 fixed discrete
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
    if(sum(dist < delta0)==0){
      m0 = m0+1
      if(d.rv==1){x = c(x, new.point)}else{x = rbind(x, new.point)}
      #x = rbind(x, new.point) #add new row of design points
      p0 = rep(1/m0,m0) #generate approx design as exp r.v.
      #update f.det
      if(MLM==TRUE){
      sum.f = 0
      for(i in 1:m0){
        if(d.rv==1){Fi = Fi.func(h.func(x[i]), bvec=bvec, link)}
        else{ Fi = Fi.func(h.func(x[i, ]), bvec=bvec, link)}
        sum.f = sum.f + (1/m0)*Fi$F_x
      } #end of for loop
      f.det = det(sum.f)
      }else{
        if(d.rv==1){m.design=length(x)} else {m.design=nrow(x);}# initial number of design points
        p.factor=length(bvec)
        X.mat = matrix(0, m.design, p.factor);  # initial model matrix X
        w.vec = rep(0, m.design);     # w vector
        for(i in 1:m.design) {
          if(d.rv==1) {htemp=Xw_maineffects_self(x=x[i], b=bvec, link=link, h.func=h.func)} else {
            htemp=Xw_maineffects_self(x=x[i,], b=bvec, link=link, h.func=h.func);
          };
          X.mat[i,]=htemp$X;
          w.vec[i]=htemp$w;
        };
        f.det=det(t(X.mat * (p0*w.vec)) %*% X.mat)    # X^T W X
      }
    } #end of if
    iter = iter + 1
  } #end of while loop


  #sort the order
  if(d.rv==1){X = x[do.call(order, as.data.frame(x))]}else
  {X = x[do.call(order, as.data.frame(x)), ]}

  list(X = X, p0 = p0, f.det=f.det)
}#end of function
