## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ForLion)
library(psych)

## -----------------------------------------------------------------------------
hfunc.temp = function(y) {c(y,y[4]*y[5],1);};   # y -> h(y)=(y1,y2,y3,y4,y5,y4*y5,1)
n.factor.temp = c(0, 2, 2, 2, 2)  # 1 continuous factor with 4 discrete factors
factor.level.temp = list(c(25,45),c(-1,1),c(-1,1),c(-1,1),c(-1,1)) 
link.temp="logit"
beta.value = c(0.35,1.50,-0.2,-0.15,0.25,0.4,-7.5)       # continuous first and intercept last to fit hfunc.temp

## -----------------------------------------------------------------------------
ForLion_GLM_Optimal(n.factor=n.factor.temp, factor.level=factor.level.temp, hfunc=hfunc.temp,bvec=beta.value,link=link.temp,reltol=1e-8, rel.diff=0.03, maxit=500, random=FALSE, logscale=TRUE)

## -----------------------------------------------------------------------------
link.temp = "cumulative"
## Note: Always put continuous factors ahead of discrete factors, pay attention to the order of coefficients paring with predictors
n.factor.temp = c(0,0,0,0,0,2)  # 1 discrete factor w/ 2 levels + 5 continuous
factor.level.temp = list(c(-25,25), c(-200,200),c(-150,0),c(-100,0),c(0,16),c(-1,1))
J = 5 #num of response levels
p = 10 #num of parameters

hfunc.temp = function(y){
  if(length(y) != 6){stop("Input should have length 6");}
  model.mat = matrix(NA, nrow=5, ncol=10, byrow=TRUE)
  model.mat[5,]=0
  model.mat[1:4,1:4] = diag(4)
  model.mat[1:4, 5] =((-1)*y[6])
  model.mat[1:4, 6:10] = matrix(((-1)*y[1:5]), nrow=4, ncol=5, byrow=TRUE)
  return(model.mat)
}
hprime.temp=NULL #use numerical gradient for optim, thus could be NULL, if use analytical gradient then define hprime function
b.temp = c(-1.77994301, -0.05287782,  1.86852211, 2.76330779, -0.94437464, 0.18504420,  -0.01638597, -0.03543202, -0.07060306, 0.10347917)

## -----------------------------------------------------------------------------
ForLion_MLM_Optimal(J=J, n.factor=n.factor.temp, factor.level=factor.level.temp, hfunc=hfunc.temp, h.prime=hprime.temp, bvec=b.temp, link=link.temp, Fi.func=Fi_MLM_func, delta=1e-2, epsilon=1e-12, reltol=1e-10, rel.diff=5e-4, maxit=500, optim_grad=FALSE)

