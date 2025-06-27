#' function for calculating X=h(x) and E_w=E(nu(beta^T h(x))) given a design point x=(1,x1,...,xd)^T
#'
#' @param x            x=(x1,...,xd) -- design point/experimental setting
#' @param Integral_based TRUE or FALSE, if TRUE then we will find the integral-based EW D-optimality otherwise we will find the sample-based EW D-optimality
#' @param b_matrix     matrix of bootstrapped or simulated parameter values.
#' @param joint_Func_b prior distribution function of model parameters
#' @param Lowerbounds  vector of lower ends of ranges of prior distribution for model parameters.
#' @param Upperbounds  vector of upper ends of ranges of prior distribution for model parameters.
#' @param link         link = "logit"  -- link function, default: "logit", other links: "probit", "cloglog", "loglog", "cauchit", "log"
#' @param h.func       function h(x)=(h1(x),...,hp(x)), default (1,x1,...,xd)
#'
#' @return X=h(x)=(h1(x),...,hp(x)) -- a row for design matrix
#' @return E_w -- E(nu(b^t h(x)))
#' @return link -- link function applied
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
#' x.temp = c(2,1,3)
#' EW_Xw_maineffects_self(x=x.temp,Integral_based=TRUE,joint_Func_b=gjoint_b,
#' Lowerbounds=paras_lowerbound,Upperbounds=paras_upperbound, link=link.temp,
#' h.func=hfunc.temp)


EW_Xw_maineffects_self <- function(x,Integral_based,joint_Func_b,Lowerbounds, Upperbounds,b_matrix,link="logit", h.func=NULL) {
  if(is.null(h.func)) h.func = function(y) {c(1,y);}; # default: main-effects
  xrow = h.func(x);
  if(Integral_based==TRUE){
  integrand_w<-function(b){
    eta = sum(b*xrow);
    w = NULL ;
    if(link=="probit") w = nu_probit_self(eta);
    if(link=="cloglog") w = nu_loglog_self(eta);
    if(link=="loglog") w = nu_loglog_self(eta);
    if(link=="cauchit") w = nu_cauchit_self(eta);
    if(link=="log") w = nu_log_self(eta);
    if(is.null(w)) { link="logit"; w=nu_logit_self(eta);};
    return(w * joint_Func_b(b))
  }
  result <- cubature::hcubature(f = integrand_w,lowerLimit = Lowerbounds,upperLimit = Upperbounds, tol = 1e-4,maxEval = 1e4)
  Ew<-result$integral
}else{
  nsa=dim(b_matrix)[1] #nsa: the number of bootstrap parameters
  w = rep(0,nsa);
  for(i in 1:nsa){
  b=b_matrix[i, ]
  eta = sum(b*xrow);
  if(link=="probit") w[i] = nu_probit_self(eta);
  if(link=="cloglog") w[i] = nu_loglog_self(eta);
  if(link=="loglog") w[i] = nu_loglog_self(eta);
  if(link=="cauchit") w[i] = nu_cauchit_self(eta);
  if(link=="log") w[i] = nu_log_self(eta);
  if(link=="logit") w[i]=nu_logit_self(eta);
  }
  Ew=mean(w,na.rm = TRUE)
}
  list(X=xrow, E_w=Ew, link=link);
}
