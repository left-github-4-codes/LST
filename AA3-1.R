#AA3 for LST by Y. Zuo on 01/13/22
#will have to call functions get_initial_beta; get_weight_wi; get_objective_value defined in AA1 and AA2

library(compiler)
enableJIT(1)
get_beta_from_p_points1=function(z)
{ # input: 
  #    z: is given data matrix with n rows and p columns, 
  #       the last column represents the y coordinates. combine 1 with (p-1) vector to get p vector x
  # output:
  #    beta [determined by p points ((x_i)', y_i): y_i=(x_i)'beta]
  
  n = dim(z)[1]
  p = dim(z)[2]
  
  beta = rep(0, p); z_temp = matrix(0, nrow=p, ncol=p)
  x_temp = matrix(0, nrow=p, ncol=p)
  
  z_temp = z[sample(1:n,p),]          #sampling p rows from z
  x_temp = cbind(1, z_temp[, 1:(p-1)])#combine 1 with x_i to form p vector and take care of intercept term
              
  
  while(det(x_temp==0)) 
    {
    z_temp = z[sample(1:n,p),]          #sampling p rows from z
    x_temp = cbind(1, z_temp[, 1:(p-1)])#combine 1 with x_i to form p vector,taking care of intercept term
    }
  
  y_temp = z_temp[,p]                   #corresponding p y_i in the equation y_i=(x_i)'beta
  
  beta = c(solve(x_temp)%*%y_temp)      #convert to a row vector
  
  return(beta)
}
get_beta_from_p_points=cmpfun(get_beta_from_p_points1)   #speed up when repeatedly call this function

AA3_main=function(z, alpha, c, N, epsilon)
{ #z is a data matrix n by p
  # alpha is used in LST
  # c is the parameter used in get_initial_beta
  # N is the total number allow for the loop
  # epsilon is the cut_off value to break loop
  
  n = dim(z)[1]
  p = dim(z)[2]
  
  a=(300)*(p-1); b=choose(n,p); N=min(a,b,N)
  wi=rep(0, p); k=0
  Q_old=10^6
  beta_old= get_initial_beta(z, c) #or use rep(0, p) instead #initial beta

  beta_new=get_beta_from_p_points(z)
  norm_diff=norm(beta_new-beta_old, type="2")
  
  while ((k<=N)& (norm_diff>=epsilon))
  { 
    beta_old=beta_new
    beta_new=get_beta_from_p_points(z)
    
    wi=get_weight_wi(z, beta_new, alpha)
    Q_new=get_objective_value(z, wi, beta_new)
    
    if(Q_new< Q_old)
    {
      Q_old=Q_new
      beta_old=beta_new
    }
    else if(Q_new==Q_old){break} #this part is modified on 01/18/22 to add this break, old one no if else
         else{
              z_temp=z[as.logical(wi),]  #pick I(beta_new) rows of z, see paper section 2.3
    
              fit1=lm(z_temp[,p]~z_temp[,1:(p-1)])
              beta_new=as.numeric(fit1$coefficients)
    
              norm_diff=norm(beta_new-beta_old, type="2")
              wi=get_weight_wi(z, beta_new, alpha)
              Q_old=get_objective_value(z, wi, beta_new)
              k=k+1
            }
  } # end while loop 

  return(beta_new)
} #end the function

# Example to call AA3_main
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
# #lst_beta=AA1_main(m1,1,0, 500, 0.001)
# lst_beta_3=AA2_main(m1,1, 0, 1000, 0.001)
# #abline(lst_beta,lty=3, col=4, lwd=1)
# abline(lst_beta_3,lty=4, col=5, lwd=1)