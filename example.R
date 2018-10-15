
rm(list=ls())

library(survAUC)


source('StrInteraction_fun.R')


m=100
p=5000
q=5

G_data_class='SNP'
model='linear'
rho=0.3
corr_str='AR'
cc_rate=20
maf=0.05


if (model=='linear'){
  n=200
} else {
  n=250
}

sigma_GE=corr_setting(p,q,corr_str,rho)





print('generate dataset')
  
data<-simulated_data(n,q,p,G_data_class,model,sigma_GE$sigmaG,sigma_GE$sigmaE,maf,cc_rate)
  
  
  
alpha_true=data$para_true$alpha
b_true=data$para_true$beta
b_vector=matrix(b_true,(q+1)*p,1)
  
G=data$G
E=data$E
y=data$Y

  
  
  
if (model=='AFT'){
  y_s=y[,1]
  delta=y[,2]
  y=y[,1]
  kmweight=kmw(y,delta)
  } else {
    kmweight=NULL
  }

gamma1=gamma2=3
lambda1=lambda2=0.12
lambda3=lambda4=0.25


result<-Spline_MCP_Hier(G,E,y,lambda4,lambda3,lambda2,lambda1,gamma1,gamma2,penalty='spline',kmweight)

alpha_estimate=result$alpha
beta_estimate=result$beta
intercept_estimate=result$intercept
beta_temp=result$beta

par(mfrow=c(2,3))

A=GetFPTP(b_vector,matrix(beta_temp,p*(q+1),1))


plot(1:100,beta_temp[1,1:100],ylim=c(-0.1,max(b_true[1,])+0.2),main=paste("TP=",A$TP," FP=",A$FP))
lines(1:100,b_true[1,1:100],'lty'=2,col='red')
lines(1:100,beta_temp[1,1:100],'lty'=1,col='blue')

plot(1:100,beta_temp[2,1:100],ylim=c(-0.1,max(b_true[2,])+0.2))
lines(1:100,b_true[2,1:100],'lty'=2,col='red')
lines(1:100,beta_temp[2,1:100],'lty'=1,col='blue')

plot(1:100,beta_temp[3,1:100],ylim=c(-0.1,max(b_true[3,])+0.4))
lines(1:100,b_true[3,1:100],'lty'=2,col='red')
lines(1:100,beta_temp[3,1:100],'lty'=1,col='blue')

plot(1:100,beta_temp[4,1:100],ylim=c(-0.1,max(b_true[4,])+0.2))
lines(1:100,b_true[4,1:100],'lty'=2,col='red')
lines(1:100,beta_temp[4,1:100],'lty'=1,col='blue')

plot(1:100,beta_temp[5,1:100],ylim=c(-0.1,0.1))
lines(1:100,b_true[5,1:100],'lty'=2,col='red')
lines(1:100,beta_temp[5,1:100],'lty'=1,col='blue')

plot(1:100,beta_temp[6,1:100],ylim=c(-0.1,0.1))
lines(1:100,b_true[6,1:100],'lty'=2,col='red')
lines(1:100,beta_temp[6,1:100],'lty'=1,col='blue')
GetFPTP(b_vector,matrix(beta_temp,p*(q+1),1))

            
          