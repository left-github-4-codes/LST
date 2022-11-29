#code for table 1 for LST by Y.Zuo on 01/14/22, modified on 11/30/22 by Y.Zuo for new table for AA1
#initilization
library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)

R=1000
n=200
p=10
alpha=1
c=0
cut_off=10^{-3}
N_ls=300
delta=1
#gamma=10^2
beta_aa1=beta_aa2=beta_aa3=matrix(0, nrow=R, ncol=p)
t_aa1=t_aa2=t_aa3=0
EMSE_aa1=EMSE_aa2=EMSE_aa3=0

#epsilon=0.05; n1=floor(epsilon*n) #consuder the case of contamination with diffeent rates of contamination

for (i in 1:R)
{# m1=rmvnorm(n, mean=(rep(0, p)),sigma=matrix(c(1, 0.88, 0.88, 1), ncol=2, byrow=T))
  z=rmvnorm(n, mean=(rep(0, p)), sigma=diag(rep(1,p)))
 # if(n1>=1)
 #  {
 #    m2=rmvnorm(n1, mean=c(7,-2),sigma=diag(rep(0.1, p)))
 #    m1[sample(1:n,n1),]<-m2 
 #  }
  
   t1=Sys.time()  
   beta_aa1[i,]=AA1_main_new_lst(z,alpha, delta, N_ls)
   t2=Sys.time()-t1
   t_aa1=t_aa1+t2
  
   # t1=Sys.time()  
   # beta_aa2[i,]=AA2_main(z,alpha, c,gamma, cut_off )
   # t2=Sys.time()-t1
   # t_aa2=t_aa2+t2
  # 
   t1=Sys.time()  
   beta_aa3[i,]=AA3_main(z,alpha, c,N, cut_off )
   t2=Sys.time()-t1
   t_aa3=t_aa3+t2
}

 EMSE_aa1=sum(beta_aa1*beta_aa1)/R
# EMSE_aa2=sum(beta_aa2*beta_aa2)/R
 EMSE_aa3=sum(beta_aa3*beta_aa3)/R

print(list("R-n-p-c-N_ls-alpha-delta-epsilon-cut_off",c(R, n,p,c, N_ls,alpha,delta,epsilon,cut_off) ))
print(list("EMSE_aa1-EMSE-aa3", c(EMSE_aa1, EMSE_aa3)))
print(list("t_aa1-t_aa3", c(t_aa1,t_aa3)))