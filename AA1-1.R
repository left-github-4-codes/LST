#AA1 for LST by Y. Zuo on 01/06/22
library(mvtnorm)
library(robustbase)

get_initial_beta=function(z,c)
{ # input: Z is a n by p matrix, first (p-1) columns for x coordinates, pth for the y coordinates
  # output: a vector beta which might be ls or lts estimator, or any pre-determined one
  # c=0, beta=rep(0,p); c=1, beta is given by LS; c=2, beta is given by LTS 
  
  p=dim(z)[2]
  
  if (c==1)
   { # get a LTS line
    fit1<-ltsReg(z[,1:(p-1)], z[,p])
    beta=as.numeric(fit1$coefficients) 
   }  
  if (c==2)
   { # get a LS line
    fit2<-lm(z[,2]~z[,1]) #p is assume to be 2
    beta=as.numeric(fit2$coefficients)
   
  }
  if (c==0)  
   { # get a pre-determined beta such as rep(0,p) or rep(1,p), etc.
    beta=rep(2,p)
   }
  return(beta)
}
#
get_weight_wi=function(z, beta, alpha)
{ # Input: Z (n by p) matrix, beta (p by 1) vector, alpha a constant
  # Output: weight w_i=indicator(|r_i-m(r_i)|/sigma(r_i)\leq alpha)  
  
  n=dim(z)[1]; p=dim(z)[2]
  Wn=rep(0, n); rn=rep(0, n); sd=rep(1, n);index=rep(0,n)
  rn=z[,p]-as.matrix(cbind(matrix(1, nrow=n), z[, 1:(p-1)]))%*%matrix(beta, nrow=p, ncol=1) 
  #residuals vector (n by 1)
  med=median(rn)
  madd=mad(rn)
  sd=abs(rn-matrix(med, nrow=n, ncol=1))/madd #generalized deviations from median 
  index=which(sd<=alpha) # which returns indeice meeting the condition
  #create a weight vector based on index
  diagg=rep(0, n) # 1 by n vector with 0 for all elements
  diagg[index]<-1 # change entries to one based on the which returned indices
  return(diagg)   # a weight vector (1 by n) with zero or one entries
}
interative_beta=function(z, wi)
{ # input: data set z (n by p) matrix, and a 1 by n vector wi given by get_weight_wi 
  # output: beta_new=(Xnn Wnn^k Xnn')^{-1}(Xnn Wnn^k Ynn)

  n=dim(z)[1]; p=dim(z)[2] 
  xnn=as.matrix(cbind(matrix(1, nrow=n), z[,1:(p-1)]))  #X matrix n by p with 1st column 1  
  
  ynn=as.matrix(z[, p])                                 # Y vector
  wnn=diag(wi)                                          # Weight matrix
  
  m=t(xnn)%*%wnn%*%xnn
  if (det(m)==0) {beta_new=ginv(m)%*%t(xnn)%*%wnn%*%ynn}
  else
  {beta_new=solve(t(xnn)%*%wnn%*%xnn)%*%t(xnn)%*%wnn%*%ynn} # weighted LS solution
  return(c(beta_new))
}
#main part
#Parameters in the function that need to be tuned every time one runs the program
AA1_main=function(z, alpha, choice_value_4_initial_beta,replication_number, cut_off)
{ # z is a given data matrix n by p, replication_number is the number user provide
  # for the total iteration allowed for the AA1 procedure (usually 200 is enough)
  # choice_value_4_initial_beta usually have three choices, 0, 1, 2, corresponding to
  # a pre-specified initial beta^0, a LS beta, or a beta from LTS
  # cut_off is a value to stop the loop of iteration
  # alpha is the number in the defintion of LST, default is one, could be 2, 3, etc.
  
  n=dim(z)[1]
  p=dim(z)[2]
  diff=beta_old=beta_new=rep(0,p)
  
  N=replication_number          #200 is more than enough in the most cases
  c=choice_value_4_initial_beta #c=0, a pre-specified beta, eg rep(0, p), rep(1,p)
                                #c=1, a LS estimator c=2, the LTS estimator
  
  alpha=1                       #default value is one, but could be 2, or 3, etc.
                                #larger one is preferred for data containing outliers
  
  threshold= cut_off #default is 10^{-3}             # stopping criterion
  beta_old=get_initial_beta(z, c)                    # initial beta^0
  wi=get_weight_wi(z, beta_old, alpha)               # weight wi is a 1 by n row vector

  for (i in (1:N))
  {    
      
    beta_new=interative_beta(z, wi) 
      
    diff=beta_new-beta_old 
    norm_diff=norm(beta_new-beta_old, type="2")
     
    if (norm_diff<=threshold){break}
    else 
    { wi=get_weight_wi(z, beta_new, alpha) # obtain a new weight vector based on the beta_new
     beta_old=beta_new                     # update beta_old and ready for the next loop
    }
  }
  
return(beta_new)
}

# Example to call AA1_main
# 
# x=c(5, 5.5, 4, 3.5, 3, 2.5, -2)
# y=c(-.5, -.5, 6, 4, 2.4, 2, .5)
# m1=cbind(x,y)
# par(mfrow=c(1,1))
# plot(m1, pch=16, xlim=c(-3,8), ylim=c(-2, 7))
# abline(0,0, lty=1, col=1, lwd=1)
# abline(0, 1, lty=2, col=2, lwd=1)
# p=dim(m1)[2]
# fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
# beta_ltsReg=as.numeric(fit1$coefficients) 
# abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=3, lty=3, lwd=1)
# lst_beta=AA1_main(m1,1,0, 500,0.001)
# abline(lst_beta,lty=3, col=4, lwd=1)