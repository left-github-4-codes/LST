# code for figure in the example 6.2 for the highly correlated contaminated normal data points
library(mvtnorm)
library(robustbase)

c=0
n=80; p=2#
epsilon=0.2; n1=floor(epsilon*n)
m1=rmvnorm(n, mean=(rep(0, p)),sigma=matrix(c(1, 0.9, 0.9, 1), ncol=2, byrow=T))
m11=m1
par(mfrow=c(1,1));  dev.off()
par(mfrow=c(1,2))  

if(n1>=1)
{
  m2=rmvnorm(n1, mean=c(7,-2),sigma=diag(rep(0.1, p)))
  m1[sample(1:n,n1),]<-m2 
}
m12=m1
#plot(m1[,1], m1[,2], xlim=c(-5, 8), ylim=c(-4, 3), xlab="x-axis", ylab="y-axis", pch=20)
plot(m1[,1], m1[,2], xlim=c(-5, 8), ylim=c(-3, 3), xlab="x-axis", ylab="y-axis", pch=20)
fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
beta_ltsReg=as.numeric(fit1$coefficients) 
abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=1, lty=1, lwd=1)
lst_beta=AA1_main(m1, 3, 0, 200,  0.001)
lst_beta_1=AA2_main(m1, 3,0, 1000, 0.001)
abline(lst_beta,lty=2, col=2, lwd=1)
abline(lst_beta_1,lty=3, col=3, lwd=1)
text(4, -1.5, expression ('LTS'))
text(4, -.5, expression('AA2'))
text(3, 2, expression('AA1'))
legend(-4, 3, legend=c("LTS", "AA1", "AA2"),    
       col=1:3, lty=1:3, cex=0.5,
       title="Line types", text.font=2, bg='lightblue')