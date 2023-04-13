####################################
##   Weibull Nadarajah-Haghighi   ##
####################################
# defining the parameter set
b = 5.7
a = 2.7
alpha_orig = alpha = 0.5
lambda_orig = lambda = 1
# generating the random sample
n=100
set.seed(100)
u = runif(n)
x = ((1+log(1+(-log(1-u)/a)^(1/b)))^(1/alpha)-1)/lambda
# WNH log-likelihood
Loglik_WNH<- function(x){
  (n*log(a*b*alpha*lambda)
   +(alpha-1)*sum(log(1+lambda*x))
   +(b-1)*sum(log(1-exp(1-(1+lambda*x)^alpha)))
   -a*sum((exp((1+lambda*x)^alpha-1)-1)^b)-b*n+b*sum((1+lambda*x)^alpha))
}
# Baseline log-likelihood
Loglik_NH<- function(x){
  sum(log(alpha*lambda*(1+lambda*x)^(alpha-1) * exp(1-(1+lambda*x)^alpha)))
}

# varying the parameters from the baseline
alphai <- seq(from = alpha_orig-.4,to = alpha_orig+.6,length.out = 500)
lambdai <- seq(from = lambda_orig-1.6,to = lambda_orig+10,length.out = 500)

# computing the log-likelihood of both models
Loglik1<-Loglik2<-Loglik3<-Loglik4<-c()
for(i in 1:length(alphai)){
  alpha<-alphai[i]
  Loglik1[i]<-Loglik_WNH(x)
  Loglik2[i]<-Loglik_NH(x)
  alpha<-alpha_orig
  lambda<-lambdai[i]
  Loglik3[i]<-Loglik_WNH(x)
  Loglik4[i]<-Loglik_NH(x)
}

###############
##   plots   ##
###############
# width and height of the plots
h11 <- w1 <- 9 

setEPS()
# Plot of the loglik varying alpha
postscript(file = "loglik_WNH_a.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
plot(alphai,Loglik1, lty=1, cex.lab=2.5, type = "l", xlab = "x", 
     ylab = expression("Log-likelihood"), cex.axis=2, 
     ylim=c(-250,1.5),xlim=c(.2,.65), col =1, lwd = 2.0)
lines(alphai,Loglik2, lty=2, type = "l", col = 2, lwd = 2.0)
legend("topright", c(expression("WNH"),
                     c(expression("NH")
                       )),
       col = c(1,2), lty= c(1,2), lwd = c(3,3), 
       bty="n", cex = 2.5)
dev.off()

# Plot of the loglik varying lambda
postscript(file = "loglik_WNH_b.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
plot(lambdai,Loglik3, lty=1, cex.lab=2.5, type = "l", xlab = "x", 
     ylab = expression("Log-likelihood"), cex.axis=2, 
     ylim=c(-250,0),xlim=c(0,8), col =1, lwd = 2.0)
lines(lambdai,Loglik4,lty=2, type = "l", col = 2, lwd = 2.0)
legend("topright", c(expression("WNH"),
                     c(expression("NH")
                     )),
       col = c(1,2), lty= c(1,2), lwd = c(3,3), 
       bty="n", cex = 2.5)
dev.off()

