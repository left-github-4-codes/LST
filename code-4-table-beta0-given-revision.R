#code for table 7 in LST by Y. Zuo on 01/16/22
#code for example with beta_0 is given in the revision on 04/25/23
library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)
library(MASS)

RepN=1000; R=RepN
p=4   # need to be tuned by different cases to 3, 4 or 5.
n=50
alpha=3 # one or three?
c=0
cut_off=10^{-3}
N_ls=3
delta=1
beta_aa1=beta_aa2=beta_aa3=matrix(0, nrow=R, ncol=p) #aa3 is actually replace by LTS
t_aa1=t_aa2=t_aa3=0                                  #aa2 is LMS

####### three methods to compute the regression lines/hyperplanes ###########

for (i in 1:RepN)
{
  # generating data, tuning based on model I (p=3), II (p=4), and III (p=5)

  #beta_0=-2; beta_1=0.1; beta_2=1 #(p=3)
  beta_0=-2; beta_1=0.1; beta_2=1;  beta_3=5 #(p=4)
  #beta_0=50; beta_1=0.1;beta_2=-2; beta_3=15; beta_4=100 #(p=5)
  
  #beta_zero=c(beta_0, beta_1, beta_2, beta_3, beta_4) #(p=5)
  beta_zero=c(beta_0, beta_1, beta_2, beta_3)   #(p=4)
  #beta_zero=c(beta_0, beta_1, beta_2) #(p=3)
  
  #x1=matrix(rcauchy(n), nrow=n, ncol=1)
  x1=matrix(rnorm(n), nrow=n, ncol=1)
  #x1=matrix(runif(n), nrow=n, ncol=1)
  #x2=matrix(rcauchy(n), nrow=n, ncol=1)
  x2=matrix(rnorm(n), nrow=n, ncol=1)
  #x2=matrix(runif(n), nrow=n, ncol=1)
  x3=matrix(rcauchy(n), nrow=n, ncol=1)
  #x3=matrix(rnorm(n), nrow=n, ncol=1)
  #x3=matrix(runif(n), nrow=n, ncol=1)
  #x4=matrix(rcauchy(n), nrow=n, ncol=1)
  #x4=matrix(rnorm(n), nrow=n, ncol=1)
  e=matrix(rnorm(n), nrow=n, ncol=1)
  #e=matrix(rcauchy(n), nrow=n, ncol=1)
  
  #y=e+matrix(beta_0, byrow=T, nrow=n, ncol=1)+beta_1*x1+beta_2*x2+beta_3*x3+ beta_4*x4 #(p=5)
  y=e+matrix(beta_0, byrow=T,  nrow=n, ncol=1)+beta_1*x1+beta_2*x2+beta_3*x3 #(p=4) 
  #y=e+matrix(beta_0, byrow=T, nrow=n, ncol=1)+beta_1*x1+beta_2*x2 #(p=3)
  #p=length(beta_zero)
  #m1=matrix(cbind(x1, x2, x3, x4, y), byrow=T, nrow=n, ncol=p) # (p=5)
  m1=matrix(cbind(x1, x2, x3, y), byrow=T, nrow=n, ncol=p) #(p=4)
  #m1=matrix(cbind(x1, x2, y), byrow=T, nrow=n, ncol=p) #(p=3)
  
  #-----------------------------------------------------------------------
  z=m1
  
  t1=Sys.time()  
  beta_aa1[i,]=AA1_main_new_lst_v2(z, alpha, delta, N_ls) #must first run the RHS function
  t2=Sys.time()-t1
  t_aa1=t_aa1+t2
  
  t1=Sys.time()  
  fit0<-lmsreg(z[,1:(p-1)], z[,p])
  beta_aa2[i,]= as.numeric(fit0$coefficients)    #AA2_main(z,alpha, c,gamma, cut_off )
  t2=Sys.time()-t1
  t_aa2=t_aa2+t2
  #print("iteration number:"); print(i)
  #print("LMS:"); print(beta_aa2[i,])
  
  t1=Sys.time()  
  fit1<-ltsReg(z[,1:(p-1)], z[,p])
  beta_aa3[i,]=as.numeric(fit1$coefficients) 
  t2=Sys.time()-t1
  t_aa3=t_aa3+t2
}  

beta_matrix=matrix(beta_zero, byrow=T, nrow=RepN, ncol=p)
beta=beta_matrix
EMSE_aa3= sum((beta_aa3-beta)*(beta_aa3-beta))/RepN
EMSE_aa1= sum((beta_aa1-beta)*(beta_aa1-beta))/RepN
EMSE_aa2= sum((beta_aa2-beta)*(beta_aa2-beta))/RepN


print(beta_zero)
print(c(R, n,p, N_ls, alpha, delta) )
print( c(EMSE_aa1, EMSE_aa2, EMSE_aa3))
print( c(t_aa1, t_aa2, t_aa3))