#' Function to generate the Expectation of fisher information at one design point xi for multinomial logit models
#'
#' @param X_x model matrix for a specific design point x_i, X_x=h.func(xi)
#' @param bvec_matrix the matrix of the bootstrap parameter values of beta
#' @param link multinomial logit model link function name "baseline", "cumulative", "adjacent", or"continuation", default to be "continuation"
#'
#' @return F_x Fisher information matrix at x_i
#' @return EU_x U matrix for calculation the Expectation of Fisher information matrix at x_i
#' @export
#'
#' @examples
#' link.temp = "continuation"
#' xi.temp=c(80)
#' hfunc.temp = function(y){
#' matrix(data=c(1,y,y*y,0,0,0,0,0,1,y,0,0,0,0,0), nrow=3, ncol=5, byrow=TRUE)
#' }
#' X_xtemp=hfunc.temp(xi.temp)
#' bvec_bootstrap<-matrix(c(-0.2401, -1.9292, -2.7851, -1.614,-1.162,
#'                          -0.0535, -0.0274, -0.0096,-0.0291, -0.04,
#'                           0.0004,  0.0003,  0.0002,  0.0003,  0.1,
#'                          -9.2154, -9.7576, -9.6818, -8.5139, -8.56),nrow=4,byrow=TRUE)
#' EW_Fi_MLM_func(X_x=X_xtemp, bvec_matrix=bvec_bootstrap, link=link.temp)

EW_Fi_MLM_func = function(X_x, bvec_matrix, link="continuation"){
  p = dim(bvec_matrix)[2]
  bs=dim(bvec_matrix)[1] ##  bs: the number of boostrap parameters
  eta_1 =  X_x %*% bvec_matrix[1, ]
  J=length(eta_1)
  U_matrix=rep(0,J*J*bs); dim(U_matrix)=c(J,J,bs)
  for(bst in 1:bs){
    eta_xi =  X_x %*% bvec_matrix[bst, ] #eta
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
    U[1:(J-1), J] =U[J,1:(J-1)]=0
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
        if(link=="continuation"){U[s,t]=U[t,s]=0}#end of continuation
        if(link=="cumulative"){
          if(t - s == 1){U[s,t]=U[t,s]= -(sum(pi[1:s])*sum(pi[1:t]))*(1-sum(pi[1:s]))*(1-sum(pi[1:t]))*(1/pi[t])}
          if(t - s > 1){U[s,t]=U[t,s]=0}
        }#end of cumulative
        if(link=="baseline"){U[s,t]=U[t,s]=-pi[s]*pi[t]}#end of baseline
        if(link=="adjacent"){U[s,t]=U[t,s]=sum(pi[1:s])*(1-sum(pi[1:t]))}
      }
    }# end of double for loop
    U_matrix[ , ,bst]=U
  }
  #EW: mean(U)
  EWU<-matrix(data = NA, nrow= J, ncol = J)
  ffi<-rep(NA,bs)
  for (i in 1:J) {
    for (j in 1:J){
      for(e in 1:bs){
        ffi[e]<-U_matrix[i,j,e]
      }
      EWU[i,j]<-mean(ffi)
    }

  }
  F_x = t(X_x) %*% EWU %*% X_x
  list(F_x = F_x, EU_x = EWU)
} #end of EW Fi.func
