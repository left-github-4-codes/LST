#code for real data sets by Y.Zuo on 01/19/22

library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)

library(MASS)

data(coleman)  #tuning the data set as you wanted
#data(wood)
#data(aircraft)
#data(salinity)
#data(delivery)
#data(phosphor)
#data(heart)


z=coleman
#z=wood
#z=aircraft
#z=salinity
#z=delivery
#z=phosphor
#z=heart
z=data.matrix(z)
n=dim(z)[1]; p=dim(z)[2]

RepN=1000; R=RepN
alpha=3 # one or three?
c=0
cut_off=10^{-3}
N_ls=10
delta=1
#gamma=10^2
beta_aa1=beta_aa2=beta_aa3=matrix(0, nrow=R, ncol=p) #aa3 is actually replace by LTS
t_aa1=t_aa2=t_aa3=0
#i=1
for (i in 1:RepN)
{
#z=m1
t1=Sys.time()  
beta_aa1[i,]=AA1_main_new_lst_v2(z,alpha, delta, N_ls)
 # AA1_main(z,alpha, c,N, cut_off )
t2=Sys.time()-t1
t_aa1=t_aa1+t2

# t1=Sys.time()  
# beta_aa2[i,]=AA2_main(z,alpha, c,gamma, cut_off )
# t2=Sys.time()-t1
# t_aa2=t_aa2+t2

t1=Sys.time()  
fit1<-ltsReg(z[,1:(p-1)], z[,p])
beta_aa3[i,]=as.numeric(fit1$coefficients) 
t2=Sys.time()-t1
t_aa3=t_aa3+t2
}


beta_aa1_mean=colMeans(beta_aa1)
deviat_beta_aa1=beta_aa1-matrix(beta_aa1_mean, nrow=RepN, ncol=p)
#beta_aa2_mean=colMeans(beta_aa2)
#deviat_beta_aa2=beta_aa2-matrix(beta_aa2_mean, nrow=RepN, ncol=p)
beta_aa3_mean=colMeans(beta_aa3)
deviat_beta_aa3=beta_aa3-matrix(beta_aa3_mean, nrow=RepN, ncol=p)

EMSE_aa1=sum(deviat_beta_aa1*deviat_beta_aa1)/RepN
#EMSE_aa2=sum(deviat_beta_aa2*deviat_beta_aa2)/RepN
EMSE_aa3=sum(deviat_beta_aa3*deviat_beta_aa3)/RepN

print(c(R, n,p,c, N_ls,alpha, delta) )
print(c(t_aa1, t_aa3))
print(c(EMSE_aa1, EMSE_aa3))

