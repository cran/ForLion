#' function to generate a initial EW Design for multinomial logistic models
#'
#' @param k.continuous number of continuous variables
#' @param factor.level lower, upper limit of continuous variables, and discrete levels of categorical variables, continuous factors come first
#' @param xlist_fix the restricted discrete settings to be chosen, default to NULL, if NULL, will generate a discrete uniform random variables
#' @param lvec lower limit of continuous variables
#' @param uvec upper limit of continuous variables
#' @param bvec_matrix the matrix of the sampled parameter values of beta
#' @param h.func function, is used to transfer the design point to model matrix (e.g. add interaction term, add intercept)
#' @param link link function, default "continuation", other options "baseline", "adjacent" and "cumulative"
#' @param EW_Fi.func function, is used to calculate the Expectation of Fisher information for a design point - default to be EW_Fi_MLM_func() in the package
#' @param delta tuning parameter, the distance threshold, || x_i(0) - x_j(0) || >= delta
#' @param epsilon determining f.det > 0 numerically, f.det <= epsilon will be considered as f.det <= 0
#' @param maxit maximum number of iterations
#'
#' @return X        matrix of initial design point
#' @return p0       initial random approximate allocation
#' @return f.det    the determinant of the expected Fisher information matrix for the initial design.
#' @export
#'
#' @examples
#' k.continuous.temp=1
#' link.temp = "continuation"
#' n.factor.temp = c(0)
#' factor.level.temp = list(c(80,200))
#' hfunc.temp = function(y){
#' matrix(data=c(1,y,y*y,0,0,0,0,0,1,y,0,0,0,0,0), nrow=3, ncol=5, byrow=TRUE)
#' }
#' lvec.temp = 80
#' uvec.temp = 200
#' bvec_bootstrap<-matrix(c(-0.2401, -1.9292, -2.7851, -1.614,-1.162,
#'                          -0.0535, -0.0274, -0.0096,-0.0291, -0.04,
#'                           0.0004,  0.0003,  0.0002,  0.0003,  0.1,
#'                          -9.2154, -9.7576, -9.6818, -8.5139, -8.56),nrow=4,byrow=TRUE)
#' EW_design_initial_MLM(k.continuous=k.continuous.temp, factor.level=n.factor.temp,xlist_fix=NULL,
#' lvec=lvec.temp,uvec=uvec.temp, bvec_matrix=bvec_bootstrap, h.func=hfunc.temp, link=link.temp)




EW_design_initial_MLM<- function(k.continuous, factor.level, xlist_fix=NULL, lvec, uvec, bvec_matrix, h.func,link="continuation", EW_Fi.func=EW_Fi_MLM_func, delta=1e-6, epsilon=1e-12, maxit=1000){
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
      sum.f = 0
      for(i in 1:m0){
        if(d.rv==1){Fi = EW_Fi.func(h.func(x[i]), bvec_matrix=bvec_matrix, link=link)}
        else{ Fi = EW_Fi.func(h.func(x[i, ]), bvec_matrix=bvec_matrix, link=link)}
        sum.f = sum.f + p0[i]*Fi$F_x
      } #end of for loop
      f.det = det(sum.f)
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
