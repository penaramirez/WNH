####################################
##   Weibull Nadarajah-Haghighi   ##
####################################
# quantile function
Q_WNH<- function(param,u){
  a=param[1]
  b=param[2]
  alpha=param[3]
  lambda=param[4]
  ((1+log(1+(-log(1-u)/a)^(1/b)))^(1/alpha)-1)/lambda
}
# moments
mh_WNH<-function(x,h){
  x^h*a*b*lambda*alpha*(1+lambda*x)^(alpha-1)*(1-exp(1-(1+lambda*x)^alpha))^(b-1)*
    exp(-b*(1-(1+lambda*x)^alpha)-a*(exp((1+lambda*x)^alpha-1)-1)^b)
}
# Bowley skewness
B_skew <- function(a,b,alpha,lambda){
  param=c(a,b,alpha,lambda)
  (Q_WNH(param,3/4)-2*Q_WNH(param,1/2)+Q_WNH(param,1/4))/
    (Q_WNH(param,3/4)-Q_WNH(param,1/4))
}
# Moors kurtosis
M_kurt <-function(a,b,alpha,lambda){
  param=c(a,b,alpha,lambda)
  (Q_WNH(param,7/8)-Q_WNH(param,5/8)+Q_WNH(param,3/8)-Q_WNH(param,1/8))/
    (Q_WNH(param,3/4)-Q_WNH(param,1/4))
}

###########
# Table 1 #
###########
param<-matrix(rep(NA,4*16),ncol=4) 
param[1,]<-c(0.8, 2.5, 2, .5)
param[2,]<-c(1.5, 2.5, 2, .5)
param[3,]<-c(2, 2.5, 2, .5)
param[4,]<-c(3, 2.5, 2, .5)
param[5,]<-c(1.5, 1.5, 2, .4)
param[6,]<-c(1.5, 2, 2, .4)
param[7,]<-c(1.5, 2.5, 2, .4)
param[8,]<-c(1.5, 3, 2, .4)
param[9,]<-c(1.4, 0.8, 0.8, .5)
param[10,]<-c(1.4, 0.8, 1.3, .5)
param[11,]<-c(1.4, 0.8, 2, .5)
param[12,]<-c(1.4, 0.8, 3, .5)
param[13,]<-c(2, 1.5, 1.8, .5)
param[14,]<-c(2, 1.5, 1.8, .8)
param[15,]<-c(2, 1.5, 1.8, 1.2)
param[16,]<-c(2, 1.5, 1.8, 1.5)

n<-dim(param)[1]
momen<-matrix(NA,nrow=n,ncol=10) # matriz que guardar?a todos os resultados 
colnames(momen)<-c("a","b","alpha","lambda","E(X)", "E(X2)","E(X3)","E(X4)","B","KM") 

mh=matrix(NA,6)
for(i in 1:n){
  a=param[i,1]
  b=param[i,2]
  alpha=param[i,3]
  lambda=param[i,4]
  for(j in 1:4){
    mh_j<-integrate(mh_WNH,0,Inf,h=j)
    mh[j]<-mh_j$value
  }
  # Bowley skewness
  mh[5]<- B_skew(a,b,alpha,lambda)
  # Moors kurtosis
  mh[6]<- M_kurt(a,b,alpha,lambda)
  momen[i,]<-c(a,b,alpha,lambda,mh)
}

print(momen, quote = FALSE, digits = 4)
