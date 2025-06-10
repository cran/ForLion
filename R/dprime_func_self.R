#' function to calculate du/dx in the gradient of d(x, Xi), will be used in ForLion_MLM_func() function, details see Appendix C in Huang, Li, Mandal, Yang (2024)
#' @param xi a vector of design point
#' @param bvec parameter of the multinomial logistic regression model
#' @param h.func function, is used to transfer xi to model matrix (e.g. add interaction term, add intercept)
#' @param h.prime function, is used to find dX/dx
#' @param inv.F.mat inverse of F_Xi matrix, inverse of fisher information of current design w/o new point
#' @param Ux U_x matrix in the algorithm, get from Fi_MLM_func() function
#' @param link multinomial link function, default is"continuation", other choices "baseline", "cumulative", and "adjacent"
#' @param k.continuous number of continuous factors
#'
#' @return dU/dx in the gradient of sensitivity function d(x, Xi)
#' @export
#'


dprime_func_self <- function(xi, bvec, h.func, h.prime, inv.F.mat, Ux, link="continuation",k.continuous){
  #step1: calculate t(C) inv(D) L
  X_x = h.func(xi)
  p = length(bvec)
  eta_xi =  X_x %*% bvec #eta
  J = length(eta_xi)
  pi = rep(NA, J) #calculate pi
  gamma = rep(NA, J-1)
  product_element = 1/(exp(eta_xi)+1)
  mat.temp = matrix(1, nrow=(J-1), ncol=(J-1))
  mat.temp[lower.tri(mat.temp)]=0
  Dj_element = mat.temp%*%eta_xi[1:(J-1)]
  Dj = sum(exp(Dj_element))+1
  #calculate pi
  for(i in 1:(J-1)){
    if(link=="continuation"){pi[i]=exp(eta_xi[i])*prod(product_element[1:i])} #end of continuation
    if(link=="cumulative"){
      if(i==1){pi[i]=exp(eta_xi[i])/(1+exp(eta_xi[i]))}
      else{
        pi[i]=(exp(eta_xi[i])/(1+exp(eta_xi[i]))) - (exp(eta_xi[i-1])/(1+exp(eta_xi[i-1])))
      }
    }#end of cumulative
    if(link=="baseline"){pi[i]=exp(eta_xi[i])/(sum(eta_xi[1:(J-1)])+1)}#end of baseline
    if(link=="adjacent"){pi[i]=(exp(Dj_element[i]))/(Dj)}#end of adjacent
  }#end of for loop
  if(link=="continuation"){pi[J]=prod(product_element[1:(J-1)])}
  if(link=="cumulative"){pi[J]=1/(1+exp(eta_xi[J-1]))}
  if(link=="baseline"){pi[J]=1/(sum(exp(eta_xi))+1)}
  if(link=="adjacent"){pi[J]=1/Dj}
  gamma = cumsum(pi) #calculate gamma
  E = diag(J)
  #calculate c1 to cJ in CDL inverse
  CDL = matrix(nrow=J, ncol=J)
  if(link=="continuation"){
    for(i in 1:J){
      if(i==1){CDL[,i]=pi[i]*c(1-gamma[i], -pi[2:J])}
      else if(i==J){CDL[,J]=pi}
      else{CDL[,i]=(pi[i]/(1-gamma[i-1]))*c(rep(0, (i-1)), 1-gamma[i], -pi[(i+1):J])}
    }
  }#end of continuation
  if(link=="cumulative"){
    for(i in 1:(J-1)){
      CDL[,i] = gamma[i]*(1-gamma[i])*(E[,i]-E[,(i+1)])
    }
    CDL[,J] = pi
  }#end of cumulative
  if(link=="baseline"){
    for(i in 1:(J-1)){
      CDL[,i]=pi[i]*(E[i]-pi)
    }
    CDL[,J]=pi
  }#end of baseline
  if(link=="adjacent"){
    for(i in 1:(J-1)){
      CDL[,i]=c((1-gamma[i])*pi[1:i], -gamma[i]*pi[(i+1):J])
    }
    CDL[,J]=pi
  }#end of adjacent

  #step2: matrix du/dpi
  du_dpi = rep(NA, J*J*J)
  dim(du_dpi) = c(J, J, J)
  du_dpi[,J,]=0
  #du_dpi[,J,] = matrix(0,nrow=J, ncol=J)
  #define du_dpi

  for(s in 1:(J-1)){
    if(link=="continuation"){
      if(s==1){du_dpi[s,s,] = c((1-gamma[1]^2), rep(pi[1]^2, J-1))}
      else{
        du_dpi[s,s,]=c(rep(0,s-1), (1-gamma[s])^2/(1-gamma[s-1])^2, rep(((pi[s]^2)/(1-gamma[s-1])^2), J-s))
      }}#end of continuation
    if(link=="cumulative"){
      du_dpi[s,s,]=Ux[s,s]*(c((2/gamma[s])*rep(1,s), (2/(1-gamma[s]))*rep(1,J-s))-(pi[s+1]/(pi[s]*(pi[s]+pi[s+1])))*E[,s]-(pi[s]/(pi[s+1]*(pi[s]+pi[s+1])))*E[,(s+1)] )
    }#end of cumulative
    if(link=="baseline"){
      du_dpi[s,s,]=c(rep(pi[s], s-1), (1-pi[s]), rep(pi[s], J-s))
    }#end of baseline
    if(link=="adjacent"){
      du_dpi[s,s,]=c(rep((1-gamma[s]), s), rep(gamma[s], J-s))
    }#end of adjacent
  }
  # upper triangle + symmetric
  for(t in 1:(J-1)){
    for(s in 1:(t-1)){
      if(link=="continuation"){du_dpi[s,t,]=c(rep(0,J))}#end of continuation
      if(link=="cumulative"){
        if((t-s)==1){
          du_dpi[s,t,]=c((gamma[s]-1)*(1-gamma[t])*(1+2*(gamma[s]/pi[t]))*rep(1,s),gamma[s]*(gamma[t]-1)*(1-((gamma[s]*(1-gamma[t]))/(pi[t]^2))), -gamma[s]*gamma[t]*(1+((2*(1-gamma[t]))/(pi[t])))*rep(1,J-s-1))
        }else{
          du_dpi[s,t,]=c(rep(0,J))
        }
      }#end of cumulative
      if(link=="baseline"){
        du_dpi[s,t,]=c(rep(0, s-1), -pi[t], rep(0, t-s-1), -pi[s], rep(0, J-t))
      }#end of baseline
      if(link=="adjacent"){
        du_dpi[s,t,]=c(rep((1-gamma[t]),s), rep(0, t-s), rep(gamma[s], J-t))
      }#end of adjacent
    }
  }
  # lower triangle
  for(t in 1:J){
    for(s in 1:(t-1)){
      du_dpi[t,s,]=du_dpi[s,t,]
    }
  }

  #step3: calculate dU/dx for each xi (only dX_dx different)
  dprime = rep(NA, k.continuous) #dprime is a vector of length continuous variable
  for(i in 1:k.continuous){
    dU_dx = matrix(0, nrow=J, ncol=J)
    dX_dx = h.prime(xi)[[i]]
    for(s in 1:J){
      for(t in 1:J){
        dU_dx[s,t]=du_dpi[s,t,] %*% CDL %*% dX_dx %*% bvec
      }
    }
    #step4: combine results to get dprime for each xi
    braket = t(dX_dx) %*% Ux %*% X_x + t(X_x) %*% dU_dx%*%X_x+t(X_x)%*%Ux%*%dX_dx
    dprime[i] = psych::tr(inv.F.mat %*% braket)
  }
  dprime<-as.vector(dprime)
}
