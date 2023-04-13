####################################
##   Weibull Nadarajah-Haghighi   ##
####################################
# linear combination
setwd("~/Insync/rrg1@de.ufpe.br/Google Drive - Shared with me/UFSM/Papers 2017/Weibull NH/RCE/figure")
a=.5; b=.5
k=15
y=rep(NA,k)
for(n in 1:k){
  result=matrix(NA,n+1,n+1)
  for(i in 0:n){
    for(j in 0:n){
      result[i+1,j+1]=((-1)^i*b*a^(i+1)*gamma(b*(i+1)+j+1))/
        (factorial(i)*factorial(j)*(b*(i+1)+j)*gamma(b*(i+1)+1))
    }
  }
  y[n]=(sum(result))
}

#############
# pdf plots #
#############
setEPS()
postscript(file = "WNH_coef.eps",width=9,height=9,family = "Times")
par(mar=c(5,6,4,1)+.1)
func<-function(x){x^0}#{return(1)}
plot(1:k,y, lty=6, type = "l",lwd = 2.0, ylab =expression("S"),
     xlab = expression("n"), cex.lab=2.5, cex.axis=2)
curve(func,from=0, to=k, add = T, lty=2, type = "l",
      col =2, lwd = 2.0)
dev.off()
