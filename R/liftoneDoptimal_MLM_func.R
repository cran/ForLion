#' function of liftone for multinomial logit model
#' @param m number of design points
#' @param p number of parameters in the multinomial logit model
#' @param Xi model matrix
#' @param J number of response levels in the multinomial logit model
#' @param thetavec model parameter
#' @param link multinomial logit model link function name "baseline", "cumulative", "adjacent", or"continuation", default to be "continuation"
#' @param reltol relative tolerance for convergence, default to 1e-5
#' @param maxit the number of maximum iteration, default to 500
#' @param p00 specified initial approximate allocation, default to NULL, if NULL, will generate a random initial approximate allocation
#' @param random TRUE or FALSE, if TRUE then the function will run with additional "nram" number of initial allocation p00, default to be TRUE
#' @param nram when random == TRUE, the function will generate nram number of initial points, default is 3
#'
#' @return p                reported D-optimal approximate allocation
#' @return p0               the initial approximate allocation that derived the reported D-optimal design
#' @return Maximum          the maximum of the determinant of the Fisher information matrix
#' @return Convergence      TRUE or FALSE, whether the algorithm converges
#' @return itmax            maximum iterations
#' @export
#'
#' @examples
#' m=5
#' p=10
#' J=5
#' factor_x = matrix(c(-1,-25,199.96,-150,-100,16,1,23.14,196.35,0,-100,
#' 16,1,-24.99,199.99,-150,0,16,-1,25,-200,0,0,16,-1,-25,-200,-150,0,16),ncol=6,byrow=TRUE)
#' Xi=rep(0,J*p*m); dim(Xi)=c(J,p,m)
#' hfunc.temp = function(y){
#' if(length(y) != 6){stop("Input should have length 6");}
#'  model.mat = matrix(NA, nrow=5, ncol=10, byrow=TRUE)
#'  model.mat[5,]=0
#'  model.mat[1:4,1:4] = diag(4)
#'  model.mat[1:4, 5] =((-1)*y[6])
#'  model.mat[1:4, 6:10] = matrix(((-1)*y[1:5]), nrow=4, ncol=5, byrow=TRUE)
#'  return(model.mat)
#'  }
#' for(i in 1:m) {
#' Xi[,,i]=hfunc.temp(factor_x[i,])
#' }
#' thetavec=c(-1.77994301, -0.05287782,  1.86852211, 2.76330779, -0.94437464, 0.18504420,
#' -0.01638597, -0.03543202, -0.07060306, 0.10347917)
#' liftoneDoptimal_MLM_func(m=m,p=p,Xi=Xi,J=J,thetavec=thetavec,
#' link="cumulative",p00=rep(1/5,5), random=FALSE)



liftoneDoptimal_MLM_func <- function(m, p, Xi, J, thetavec,link = "continuation", reltol=1e-5, maxit=500, p00=NULL, random=FALSE, nram=3) {
  if(is.null(p00)){p00=stats::rexp(m); p00=p00/sum(p00);}
  Fi <- rep(0, p*p*m);  dim(Fi)=c(p,p,m)
  nFi <- rep(0, p*p*m);  dim(nFi)=c(p,p,m)

  for(i in 1:m) {
    Fi[,,i]=Fi_MLM_func(Xi[, ,i], bvec=thetavec, link=link)$F_x
    nFi[,,i]=p00[i]*Fi[,,i]
  }
  F=apply(nFi,c(1,2),sum)
  Fdet=det(F)

  Bn1 <- matrix(1, J-1, J-1)     # B^{-1}
  for(j in 2:(J-1)) Bn1[,j]=(1:(J-1))^(j-1);
  Bn1 = solve(Bn1);
  fdet <- function(p) {   # |F|=|sum_i p_i F_i|, p[1:m], need "Fi"
    atemp=p[1]*Fi[,,1];
    for(i in 2:m) atemp=atemp+p[i]*Fi[,,i];
    det(atemp);
  }
  fiz <- function(z, p, i) {   # f_i(z), need "fdet"
    p1=p*(1-z)/(1-p[i]);
    p1[i]=z;
    fdet(p1);
  }
  if(is.null(p00)) p00=p00;        # default initial point is uniform design
  maximum = fdet(p00);
  maxvec = stats::rexp(m);
  convergence = F;
  p0 = p00;
  ind = 0;

  while((ind < maxit) && ((max(maxvec)/min(maxvec))-1 > reltol)) {
    io = sample(seq(1:m));
    for(ia in 1:m) {    # run updating in random order of {1,2,...,m}
      avec <- rep(0, J);    # a0, a1, ..., a_{J-1}
      avec[1] = fiz(0, p0, io[ia]);  # a0=f_i(0)
      cvec <- rep(0, J-1);  # c1, c2, ..., c_{J-1}
      for(j in 1:(J-1)) cvec[j]=(j+1)^p*j^(J-1-p)*fiz(1/(j+1), p0, io[ia])-j^(J-1)*avec[1];
      avec[J:2]=Bn1%*%cvec;
      if(J<=5){
        ftemp <- function(z) {   # f_i(z)
          obj=(1-z)^(p-J+1)*sum(avec*z^(0:(J-1))*(1-z)^((J-1):0));
          # cat("\navec", avec, "\nz",z, "\nsum", sum(avec*z^(0:(J-1))*(1-z)^((J-1):0)),"\nobject",obj,"\n") #delete
          return(obj)
        }
        if(J==3){
          c0.temp <- avec[1]*p - avec[2]
          c1.temp <- avec[2]*p+avec[2]-2*avec[1]*p-2*avec[3]
          c2.temp <- avec[1]*p - avec[2]*p + avec[3]*p
          sol.temp = polynomial_sol_J3(c0.temp, c1.temp, c2.temp) #return the two analytical solutions
        }
        if(J==4){
          c0.temp <- avec[1]*p - avec[2]
          c1.temp <- -3*avec[1]*p + 2*avec[2] + p*avec[2] - 2*avec[3]
          c2.temp <- 3*avec[1]*p - (1+2*p)*avec[2] + (2+p)*avec[3] - 3*avec[4]
          c3.temp <- p*(-avec[1] + avec[2] - avec[3] + avec[4])
          sol.temp = polynomial_sol_J4(c0.temp, c1.temp, c2.temp, c3.temp) #return the three analytical solutions
        }
        if(J==5){
          #define the coefficients of 4th order polynomial function
          c0.temp = -avec[2]+avec[1]*p
          c1.temp = 3*avec[2] - 2*avec[3] - 4*avec[1]*p + avec[2]*p
          c2.temp = -3*avec[2] + 4*avec[3] - 3*avec[4] + 6*avec[1]*p - 3*avec[2]*p + avec[3]*p
          c3.temp = avec[2] - 2*avec[3] + 3*avec[4] - 4*avec[5] -4*avec[1]*p + 3*avec[2]*p - 2*avec[3]*p + avec[4]*p
          c4.temp =avec[1]*p - avec[2]*p + avec[3]*p - avec[4]*p + avec[5]*p
          sol.temp = polynomial_sol_J5(c0.temp, c1.temp, c2.temp, c3.temp, c4.temp) #return the four analytical solutions
        }
        #remove the complex solution, only use real solution
        sol.temp[abs(Im(sol.temp)) > 1e-6] = NA
        sol.temp[Re(sol.temp) < 1e-6] = NA
        sol.temp[Re(sol.temp) > (1-1e-6)] = NA
        sol.temp = Re(stats::na.omit(sol.temp))

        #all the four solutions are complex zstar=0 if not find the max ftemp
        zstar=0; fstar=avec[1];
        if(length(sol.temp)>0){
          for(value in sol.temp){
            ftemp.value = ftemp(value)
            if(ftemp.value > fstar){zstar=value; fstar=ftemp.value}
          }#for loop end
        } #end if
      }else{
      ftemp <- function(z) {   # -f_i(z)
        obj=-(1-z)^(p-J+1)*sum(avec*z^(0:(J-1))*(1-z)^((J-1):0));
        # cat("\navec", avec, "\nz",z, "\nsum", sum(avec*z^(0:(J-1))*(1-z)^((J-1):0)),"\nobject",obj,"\n") #delete
        return(obj)
      }
      ftemp1 <- function(z) {  # -f'_i(z)
        #   -(1-z)^(p-J)*sum(((1:(J-1))-P*z)*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+p*avec[1]*(1-z)^(p-1);
        -(1-z)^(p-J+1)*sum((1:(J-1))*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+(1-z)^(p-J)*sum((p:(p-J+1))*avec[1:J]*z^(0:(J-1))*(1-z)^((J-1):0));
      }
      temp=stats::optim(par=0.5, fn=ftemp, gr=ftemp1, method="L-BFGS-B", lower=0, upper=1, control=list(maxit=maxit, factr=1e5));
      zstar=temp$par;         # z_*
      fstar=-temp$value;
      if(fstar <= avec[1]) {zstar=0; fstar=avec[1];};
      }
      ptemp1 = p0*(1-zstar)/(1-p0[io[ia]]);
      ptemp1[io[ia]] = zstar;
      if(fstar > maximum) {maximum = fstar; p0=ptemp1;};
      maxvec[io[ia]] = maximum;
    }
    ind = ind+1;
    #cat("\nmaxit", maxit, "\nmax(maxvec)", max(maxvec), "\nmin(maxvec)", min(maxvec)) #delete
  }# end of "while"
  p00.ans = p00;
  p0.ans=p0;
  maximum.ans=maximum;
  #maximum.adj=maximum*n^p;
  #fdet.adj=Fdet*n^p;
  if((max(maxvec)/min(maxvec))-1 <= reltol) convergence=T;
  itmax=ind;
  effi=(Fdet/maximum.ans)^(1/p)
  #random initial weights
  if(random){
    for(j in 1:nram){
      p00=stats::rexp(m)
      p00=p00/sum(p00)
      Fi <- rep(0, p*p*m);  dim(Fi)=c(p,p,m)
      nFi <- rep(0, p*p*m);  dim(nFi)=c(p,p,m)

      for(i in 1:m) {
        Fi[,,i]=Fi_MLM_func(Xi[, ,i], bvec=thetavec, link=link)$F_x
        nFi[,,i]=p00[i]*Fi[,,i]
      }
      F=apply(nFi,c(1,2),sum)
      Fdet=det(F)

      Bn1 <- matrix(1, J-1, J-1)     # B^{-1}
      for(j in 2:(J-1)) Bn1[,j]=(1:(J-1))^(j-1);
      Bn1 = solve(Bn1);
      fdet <- function(p) {   # |F|=|sum_i p_i F_i|, p[1:m], need "Fi"
        atemp=p[1]*Fi[,,1];
        for(i in 2:m) atemp=atemp+p[i]*Fi[,,i];
        det(atemp);
      }
      fiz <- function(z, p, i) {   # f_i(z), need "fdet"
        p1=p*(1-z)/(1-p[i]);
        p1[i]=z;
        fdet(p1);
      }
      if(is.null(p00)) p00=p00;        # default initial point is uniform design
      maximum = fdet(p00);
      maxvec = stats::rexp(m);
      convergence = F;
      p0 = p00;
      ind = 0;

      while((ind < maxit) && ((max(maxvec)/min(maxvec))-1 > reltol)) {
        io = sample(seq(1:m));
        for(ia in 1:m) {    # run updating in random order of {1,2,...,m}
          avec <- rep(0, J);    # a0, a1, ..., a_{J-1}
          avec[1] = fiz(0, p0, io[ia]);  # a0=f_i(0)
          cvec <- rep(0, J-1);  # c1, c2, ..., c_{J-1}
          for(j in 1:(J-1)) cvec[j]=(j+1)^p*j^(J-1-p)*fiz(1/(j+1), p0, io[ia])-j^(J-1)*avec[1];
          avec[J:2]=Bn1%*%cvec;
          if(J<=5){
            ftemp <- function(z) {   # f_i(z)
              obj=(1-z)^(p-J+1)*sum(avec*z^(0:(J-1))*(1-z)^((J-1):0));
              # cat("\navec", avec, "\nz",z, "\nsum", sum(avec*z^(0:(J-1))*(1-z)^((J-1):0)),"\nobject",obj,"\n") #delete
              return(obj)
            }
            if(J==3){
              c0.temp <- avec[1]*p - avec[2]
              c1.temp <- avec[2]*p+avec[2]-2*avec[1]*p-2*avec[3]
              c2.temp <- avec[1]*p - avec[2]*p + avec[3]*p
              sol.temp = polynomial_sol_J3(c0.temp, c1.temp, c2.temp) #return the two analytical solutions
            }
            if(J==4){
              c0.temp <- avec[1]*p - avec[2]
              c1.temp <- -3*avec[1]*p + 2*avec[2] + p*avec[2] - 2*avec[3]
              c2.temp <- 3*avec[1]*p - (1+2*p)*avec[2] + (2+p)*avec[3] - 3*avec[4]
              c3.temp <- p*(-avec[1] + avec[2] - avec[3] + avec[4])
              sol.temp = polynomial_sol_J4(c0.temp, c1.temp, c2.temp, c3.temp) #return the three analytical solutions
            }
            if(J==5){
              #define the coefficients of 4th order polynomial function
              c0.temp = -avec[2]+avec[1]*p
              c1.temp = 3*avec[2] - 2*avec[3] - 4*avec[1]*p + avec[2]*p
              c2.temp = -3*avec[2] + 4*avec[3] - 3*avec[4] + 6*avec[1]*p - 3*avec[2]*p + avec[3]*p
              c3.temp = avec[2] - 2*avec[3] + 3*avec[4] - 4*avec[5] -4*avec[1]*p + 3*avec[2]*p - 2*avec[3]*p + avec[4]*p
              c4.temp =avec[1]*p - avec[2]*p + avec[3]*p - avec[4]*p + avec[5]*p
              sol.temp = polynomial_sol_J5(c0.temp, c1.temp, c2.temp, c3.temp, c4.temp) #return the four analytical solutions
            }
            #remove the complex solution, only use real solution
            sol.temp[abs(Im(sol.temp)) > 1e-6] = NA
            sol.temp[Re(sol.temp) < 1e-6] = NA
            sol.temp[Re(sol.temp) > (1-1e-6)] = NA
            sol.temp = Re(stats::na.omit(sol.temp))

            #all the four solutions are complex zstar=0 if not find the max ftemp
            zstar=0; fstar=avec[1];
            if(length(sol.temp)>0){
              for(value in sol.temp){
                ftemp.value = ftemp(value)
                if(ftemp.value > fstar){zstar=value; fstar=ftemp.value}
              }#for loop end
            } #end if
          }else{
            ftemp <- function(z) {   # -f_i(z)
              obj=-(1-z)^(p-J+1)*sum(avec*z^(0:(J-1))*(1-z)^((J-1):0));
              # cat("\navec", avec, "\nz",z, "\nsum", sum(avec*z^(0:(J-1))*(1-z)^((J-1):0)),"\nobject",obj,"\n") #delete
              return(obj)
            }
            ftemp1 <- function(z) {  # -f'_i(z)
              #   -(1-z)^(p-J)*sum(((1:(J-1))-P*z)*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+p*avec[1]*(1-z)^(p-1);
              -(1-z)^(p-J+1)*sum((1:(J-1))*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+(1-z)^(p-J)*sum((p:(p-J+1))*avec[1:J]*z^(0:(J-1))*(1-z)^((J-1):0));
            }
            temp=stats::optim(par=0.5, fn=ftemp, gr=ftemp1, method="L-BFGS-B", lower=0, upper=1, control=list(maxit=maxit, factr=1e5));
            zstar=temp$par;         # z_*
            fstar=-temp$value;
            if(fstar <= avec[1]) {zstar=0; fstar=avec[1];};
          }
          ptemp1 = p0*(1-zstar)/(1-p0[io[ia]]);
          ptemp1[io[ia]] = zstar;
          if(fstar > maximum) {maximum = fstar; p0=ptemp1;};
          maxvec[io[ia]] = maximum;
        }
        ind = ind+1;
        #cat("\nmaxit", maxit, "\nmax(maxvec)", max(maxvec), "\nmin(maxvec)", min(maxvec)) #delete
      }# end of "while"
      if(maximum > maximum.ans){
        p00.ans = p00;
        p0.ans=p0;
        maximum.ans=maximum;
        #maximum.adj=maximum*n^p;
        #fdet.adj=Fdet*n^p;
        if((max(maxvec)/min(maxvec))-1 <= reltol) convergence=T;
        itmax=ind;
        #effi=(Fdet/maximum.ans)^(1/p)
      }

    }#end of for loop of 1:nram

  }#end of if(random)

  # list(p=p0.ans, p0=p00, Maximum=maximum.ans, convergence=convergence, itmax=itmax);

  #define S3 class
  output<-list(p=p0.ans, p0=p00, Maximum=maximum.ans, convergence=convergence, itmax=itmax);
  class(output) <- "list_output"
  return(output)
}
