#' function to generate fisher information at one design point xi for multinomial logit models
#'
#' @param X_x model matrix for a specific design point x_i, X_x=h.func(xi)
#' @param bvec beta coefficients in the model
#' @param link multinomial logit model link function name "baseline", "cumulative", "adjacent", or"continuation", default to be "continuation"
#'
#' @return F_x    Fisher information matrix at x_i
#' @return U_x    U matrix for calculation of Fisher information matrix at x_i (see Corollary 3.1 in Bu, Majumdar, Yang(2020))
#' @export
#'
#' @examples
#' # Reference minimizing surface example in supplementary material
#' # Section S.3 in Huang, Li, Mandal, Yang (2024)
#' xi.temp = c(-1, -25, 199.96, -150, -100, 16)
#' hfunc.temp = function(y){
#'if(length(y) != 6){stop("Input should have length 6");}
#' model.mat = matrix(NA, nrow=5, ncol=10, byrow=TRUE)
#' model.mat[5,]=0
#' model.mat[1:4,1:4] = diag(4)
#' model.mat[1:4, 5] =((-1)*y[6])
#' model.mat[1:4, 6:10] = matrix(((-1)*y[1:5]), nrow=4, ncol=5, byrow=TRUE)
#' return(model.mat)
#' }
#' X_x.temp = hfunc.temp(xi.temp)
#' bvec.temp = c(-1.77994301, -0.05287782,  1.86852211, 2.76330779, -0.94437464,
#' 0.18504420,  -0.01638597, -0.03543202, -0.07060306, 0.10347917)
#' link.temp = "cumulative"
#' Fi_MLM_func(X_x=X_x.temp, bvec=bvec.temp, link=link.temp)
#'
#'



Fi_MLM_func = function(X_x, bvec, link="continuation"){
  p = length(bvec)
  eta_xi =  X_x %*% bvec #eta
  J = length(eta_xi)
  pi = rep(NA, J) #calculate pi
  product_element = 1/(exp(eta_xi)+1)
  #adjacent denominator
  sum_element = rep(NA, J-1)
  for(i in 1:(J-1)){
    sum_element[i] = exp(sum(eta_xi[i:(J-1)]))
  }
  for(i in 1:(J-1)){
    if(link=="continuation"){pi[i]=exp(eta_xi[i])*prod(product_element[1:i])} #end of continuation
    if(link == "cumulative"){
      if(i==1){pi[i]=exp(eta_xi[i])/(1+exp(eta_xi[i]))}else{
        pi[i]=((exp(eta_xi[i]))/(1 + exp(eta_xi[i])))-((exp(eta_xi[i-1]))/(1+exp(eta_xi[i-1])))
      }
    }#end of cumulative
    if(link=="baseline"){pi[i]=(exp(eta_xi[i]))/(sum(exp(eta_xi)[1:(J-1)])+1)}#end of baseline
    if(link=="adjacent"){pi[i]=(exp(sum(eta_xi[i:(J-1)])))/(sum(sum_element)+1)}
  }#end of for loop
  if(link=="continuation"){pi[J] = prod(product_element[1:(J-1)])}
  if(link=="cumulative"){pi[J] = 1/(1+exp(eta_xi[J-1]))}
  if(link=="baseline"){pi[J] = 1/(sum(exp(eta_xi)[1:(J-1)])+1)}
  if(link=="adjacent"){pi[J] = 1/(sum(sum_element)+1)}
  #calculate U matrix, U is symmetric matrix
  U = matrix(data = NA, nrow= J, ncol = J)
  U[1:(J-1), J] = 0
  U[J, J] = 1
  #diagonal
  if(link=="continuation"){U[1,1]= pi[1]*(1-pi[1])}
  if(link=="cumulative"){U[1,1]=(pi[1]^2)*((1-pi[1])^2)*(1/pi[1] + 1/pi[2])}
  if(link=="baseline"){U[1,1]=pi[1]*(1-pi[1])}
  if(link=="adjacent"){U[1,1]=pi[1]*(1-pi[1])}
  for(s in 2:(J-1)){
    if(link=="continuation"){U[s,s]=pi[s]*((1-sum(pi[1:s])))*((1-sum(pi[1:(s-1)]))^(-1))}
    if(link=="cumulative"){U[s,s]=(sum(pi[1:s])^2)*((1-sum(pi[1:s]))^2)*(1/pi[s] + 1/pi[s+1])}
    if(link=="baseline"){U[s,s]=pi[s]*(1-pi[s])}
    if(link=="adjacent"){U[s,s]=sum(pi[1:s])*(1-sum(pi[1:s]))}
  }#end of for loop

  #upper triangle
  for(s in 1 : (J-2)){
    for(t in (s+1):(J-1)){
      if(link=="continuation"){U[s,t]=0}#end of continuation
      if(link=="cumulative"){
        if(t - s == 1){U[s,t] = -(sum(pi[1:s])*sum(pi[1:t]))*(1-sum(pi[1:s]))*(1-sum(pi[1:t]))*(1/pi[t])}
        if(t - s > 1){U[s,t]=0}
      }#end of cumulative
      if(link=="baseline"){U[s,t]=-pi[s]*pi[t]}#end of baseline
      if(link=="adjacent"){U[s,t]=sum(pi[1:s])*(1-sum(pi[1:t]))}
    }
  }# end of double for loop

  #symmetric U codes, lower triangle = transpose lower triangle (update 11/17/2024)
  ind_temp=lower.tri(U)
  U[ind_temp]=t(U)[ind_temp]

  F_x = t(X_x) %*% U %*% X_x
  list(F_x = F_x, U_x = U)
} #end of Fi.func
