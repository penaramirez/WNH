rm(list = ls())
#############
# Scenarios #
#############
scenarios<-4
param<-matrix(rep(NA,4*scenarios),ncol=4) 
# cenarios 
param[1,]<-c(6.5, 3, .1, 1)
param[2,]<-c(3, 6.5, .1, 1)
param[3,]<-c(1.5, 6.5,.1,1)
param[4,]<-c(5.7, 2.7,.1,1)

#############
# loglik #
#############
fr <- function(par,x){
  a <- abs(par[1])
  b <- abs(par[2])
  alpha <- abs(par[3])
  l <- abs(par[4])
  -(n*log(a*b*alpha*l)
    +(alpha-1)*sum(log(1+l*x))
    +(b-1)*sum(log(1-exp(1-(1+l*x)^alpha)))
    -a*sum((exp((1+l*x)^alpha-1)-1)^b)-b*n+b*sum((1+l*x)^alpha)
  )
}

Q_WNH<- function(u){
  ((1+log(1+(-log(1-u)/a)^(1/b)))^(1/alpha)-1)/lambda
}

######################
# initial quantities #
######################
start_time <- Sys.time()
set.seed(2000)
n = 60 # sample size
nrep = 10000 # Monte Carlo replications
eta<-0.05 # 100(1-eta)% confidence level
# to save the results
estim<-EQM<-inf<-array(NA, c(dim(param)[2],nrep,scenarios))
CR<-Mu<-rb<-rmse<-matrix(nrow = dim(param)[1],ncol=dim(param)[2]) 
colnames(CR)<-colnames(Mu)<-
  colnames(rb)<-colnames(rmse)<-c("a","b","alpha","lambda")
method_simu<-"Nelder-Mead"
####################
# simulation study #
####################
for(j in 1:scenarios){
  print(j)
  a=param[j,1]
  b=param[j,2]
  alpha=param[j,3]
  lambda=param[j,4]
  par<-c(a,b,alpha,lambda)
  epsilon=1e-2
  chute=c(a+epsilon,b+epsilon,alpha+epsilon,lambda+epsilon)
  iter=0
  bug<-0
  i=1
  ### Monte Carlo loop
  while (i<=nrep){
    iter=iter+1
    u = runif(n, min = 0, max = 1)
    t_i = Q_WNH(u)
    res<-try(optim(chute, fr,  method = method_simu, hessian=T,x=t_i),
             silent = T)
    if(class(res)=="try-error" || res$conv != 0) # a classe dos objetos que cont?m o erro, 
    {
      bug<-bug+1
    }else{
      inf[,i,j]=diag(solve(-res$hessian))
      if(res$conv==0 && sum(inf[,i,j]>0)==4){
        estim[,i,j]=res$par
        EQM[,i,j]=(estim[,i,j]-par)^2
        i = i + 1
        # print(c("j=",j,"i=",i))
      } 
    }
  } # end of the Monte Carlo loop
  ### Coverage rates
  li=estim[,,j]+qnorm(eta/2)*sqrt(inf[,,j])
  ls=estim[,,j]+qnorm(1-eta/2)*sqrt(inf[,,j])
  CR[j,1]<-sum(a>li[1,] & a<ls[1,])/nrep*100
  CR[j,2]<-sum(b>li[2,] & b<ls[2,])/nrep*100
  CR[j,3]<-sum(alpha>li[3,] & alpha<ls[3,])/nrep*100
  CR[j,4]<-sum(lambda>li[4,] & lambda<ls[4,])/nrep*100
  ### Mean estimates
  Mu[j,] <- apply(estim[,,j],1,mean)
  ### Bias
  rb[j,]<- (Mu[j,]-par)/par*100
  ### RMSEs
  rmse[j,]<-sqrt(apply(EQM[,,j],1,mean))
  print(CR[j,])
}

end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)
saving<-paste0("results/",method_simu,"_n",n,"R",nrep,"lambda",lambda,".RData")
save.image(saving)

print(rb)
beepr::beep(8)
