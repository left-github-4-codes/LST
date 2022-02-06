#code for Figure one (a) and (b) By Y. Zuo on 01/22/22
####Part (a)

x=c(5, 5.5, 4, 3.5, 3, 2.5, -2)
y=c(-.5, -.5, 6, 4, 2.4, 2, .5)
m1=cbind(x,y)
par(mfrow=c(1,2))
plot(m1, pch=16, xlim=c(-3,8), ylim=c(-2, 7))
abline(0,0, lty=1, col=1, lwd=1)
abline(0, 1, lty=2, col=2, lwd=1)
text(x, y-0.5, as.character(seq(1:7)))
text(7.5, -.5, expression('L'[1]))
text(7, 6, expression('L'[2]))

p=dim(m1)[2]
plot(m1, pch=16, xlim=c(-3,8), ylim=c(-2, 7))
fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
beta_ltsReg=as.numeric(fit1$coefficients) 
abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=1, lty=1, lwd=1)
lst_beta=AA1_main(m1,1,500,0, 0.001)
abline(lst_beta,lty=2, col=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=3, lty=3, lwd=1)
text(x, y-0.5, as.character(seq(1:7)))
text(7.1, 1.7, expression ('LTS'))
text(6.7, 5.7, expression('LST'))
text(7.1, 2.5, expression('LS'))
legend(-2.5, 6.5, legend=c(#expression('L'[1]),expression('L'[2]), 
                           "LTS", "LST", "LS"),    
       col=1:3, lty=1:3, cex=0.6,
       title="Line types", text.font=3, bg='lightblue')

######Part (b)
alpha=1     #cut-off value alpha
N=500               # the total iteration number allowed
c=0                #c=0,initial beta is rep(0, p);c=1,it is by LS; c=2, it is by LTS.

n=7; p=2#
epsilon=0.3; n1=floor(epsilon*n)
m1=rmvnorm(n, mean=(rep(0, p)),sigma=matrix(c(1, 0.88, 0.88, 1), ncol=2, byrow=T))

par(mfrow=c(1,2))  
plot(m1[,1], m1[,2], xlim=c(-3,8), ylim=c(-2, 7), xlab="x-axis", ylab="y-axis", pch=20)
fit0<-ltsReg(m1[,1:(p-1)], m1[,p])
abline(fit0$coefficients[[1]], fit0$coefficients[[2]], col=1, lty=1, lwd=1)
abline(AA1_main(m1,1,500, 0, 0.001), col=2, lty=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=3, lty=3, lwd=1)
legend(-2.5, 6.5, legend=c("LTS", "LST", "LS"),    
       col=1:3, lty=1:3, cex=0.6,
       title="Line types", text.font=3, bg='lightblue')

if(n1>=1)
{
  m2=rmvnorm(n1, mean=c(7,-2),sigma=diag(rep(0.1, p)))
  m1[sample(1:n,n1),]<-m2 
}

plot(m1[,1], m1[,2], xlim=c(-3,8), ylim=c(-2, 7), xlab="x-axis", ylab="y-axis", pch=20)
fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=1, lty=1, lwd=1)
abline(AA1_main(m1,1,500, 0, 0.001), col=2, lty=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=3, lty=3, lwd=1)

text(4., -1.6, expression ('LTS'))
text(6, 5.7, expression('LST'))
text(4, 0, expression('LS'))

legend(-2.5, 6.5, legend=c("LTS", "LST", "LS"),    
  col=1:3, lty=1:3, cex=0.6,
  title="Line types", text.font=3, bg='lightblue')
