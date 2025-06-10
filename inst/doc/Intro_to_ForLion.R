## ----include=FALSE------------------------------------------------------------
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
set.seed(482)
forlion.glm=ForLion_GLM_Optimal(n.factor=n.factor.temp,factor.level=factor.level.temp, hfunc=hfunc.temp,bvec=beta.value,link=link.temp,reltol=1e-8, rel.diff=0.03, maxit=500, random=FALSE, logscale=TRUE)

forlion.glm

## -----------------------------------------------------------------------------
GLM_Exact_Design(k.continuous=1, design_x=forlion.glm$x.factor, design_p=forlion.glm$p,
det.design=forlion.glm$det,p=7, ForLion=TRUE, bvec=beta.value, rel.diff=0.5,L=0.1,N=500, hfunc=hfunc.temp, link="logit")

## -----------------------------------------------------------------------------
nrun = 200
set.seed(0713)
b_0 = runif(nrun, -8, -7)
b_1 = runif(nrun, 1, 2)
b_2 = runif(nrun, -0.3, -0.1)
b_3 = runif(nrun, -0.3, 0)
b_4 = runif(nrun, 0.1, 0.4)
b_5 = runif(nrun, 0.25, 0.45)
b_34= runif(nrun, 0.35, 0.45)
beta.matrix = cbind(b_5,b_1,b_2,b_3,b_4,b_34,b_0)

## -----------------------------------------------------------------------------
set.seed(482)
EW_ForLion_GLM_Optimal(n.factor=c(0, 2, 2, 2, 2), factor.level=list(c(25,45),c(-1,1),c(-1,1),c(-1,1),c(-1,1)), hfunc=hfunc.temp, Integral_based=FALSE, b_matrix=beta.matrix, link="logit", reltol=1e-8, rel.diff=5e-3, optim_grad=TRUE, maxit=500, random=FALSE,nram=1,logscale=TRUE)

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
set.seed(123)
ForLion_MLM_Optimal(J=J, n.factor=n.factor.temp, factor.level=factor.level.temp, hfunc=hfunc.temp, h.prime=hprime.temp, bvec=b.temp, link=link.temp, Fi.func=Fi_MLM_func, delta=1e-2, epsilon=1e-10, reltol=1e-8, rel.diff=0.5, maxit=500, optim_grad=FALSE)

## -----------------------------------------------------------------------------
nrun = 100
set.seed(0713)
b_clean = runif(nrun, -1, 0)
b_temperature = runif(nrun, 0, 0.2)
b_pressure = runif(nrun, -0.1, 0.1)
b_nitrogen = runif(nrun, -0.1, 0.1)
b_silane = runif(nrun, -0.1, 0.1)
b_time = runif(nrun, 0, 0.2)
theta1 = runif(nrun, -2, -1)
theta2 = runif(nrun, -0.5, 0.5)
theta3 = runif(nrun, 1, 2)
theta4 = runif(nrun, 2.5, 3.5)

beta.temp2 = cbind(theta1, theta2, theta3, theta4, b_clean, b_temperature, b_pressure, b_nitrogen, b_silane, b_time)

## -----------------------------------------------------------------------------
set.seed(123)
EW_forlion.MLM =EW_ForLion_MLM_Optimal(J=J, n.factor=n.factor.temp,
factor.level=factor.level.temp,hfunc=hfunc.temp, h.prime=hprime.temp, bvec_matrix=beta.temp2, link=link.temp, EW_Fi.func=EW_Fi_MLM_func, 
delta=1e-2, epsilon=1e-10, reltol=1e-8, rel.diff=0.5, maxit=500, optim_grad=FALSE)

EW_forlion.MLM 

## -----------------------------------------------------------------------------
MLM_Exact_Design(J=J, k.continuous=5,design_x=EW_forlion.MLM$x.factor,design_p=EW_forlion.MLM$p,det.design=EW_forlion.MLM$det,p=10,ForLion=FALSE,bvec_matrix=beta.temp2,rel.diff=1,L=c(0.5,0.1,0.1,0.1,1),N=1000,hfunc=hfunc.temp,link=link.temp)

