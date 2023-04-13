####################################
##   Weibull Nadarajah-Haghighi   ##
####################################
setwd("~/Insync/rrg1@de.ufpe.br/Google Drive - Shared with me/UFSM/Papers 2017/Weibull NH/RCE/figure")
# moments
mh_WNH<-function(x,h){
  x^h*a*b*lambda*alpha*(1+lambda*x)^(alpha-1)*(1-exp(1-(1+lambda*x)^alpha))^(b-1)*
    exp(-b*(1-(1+lambda*x)^alpha)-a*(exp((1+lambda*x)^alpha-1)-1)^b)
}
# Bonferroni curve
Bp<-Vectorize(function(p){
    q=((1+log(1+(-log(1-p)/a)^(1/b)))^(1/alpha)-1)/lambda
    L=integrate(mh_WNH,0,q,h=1)$value/
    (p*integrate(mh_WNH,0,Inf,h=1)$value)
    return(L)
    })
# Lorenz curve
Lp<-Vectorize(function(p){
  q=((1+log(1+(-log(1-p)/a)^(1/b)))^(1/alpha)-1)/lambda
  L=integrate(mh_WNH,0,q,h=1)$value/
    integrate(mh_WNH,0,Inf,h=1)$value
  return(L)
})

##########################
# Plots of B(p) and L(p) #
##########################
h11 <- w1 <- 9 # width and height of the plots
select<-c(1,6,12,16)
setEPS()
postscript(file = "Lp.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
a=0.8; b=2.5; alpha=2; lambda=.5
curve(Lp,from=0.00095, to=1, 
      add = FALSE, lty=1, cex.lab=2.5, type = "l", xlab = "p",
      ylab = expression("L(p)"), cex.axis=2, ylim =c(0,1.22), lwd = 2.0)
a=1.5; b=2; alpha=2; lambda=.4
curve(Lp,from=0.00095, to=1, col=(2),
        add = T, lty=2, type = "l",  lwd = 2.0)
a=1.4; b=.8; alpha=3; lambda=.5
curve(Lp,from=0.00095, to=1, col=(3),
      add = T, lty=3, type = "l",  lwd = 6.0)
a=2; b=1.5; alpha=1.8; lambda=1.5
curve(Lp,from=0.00095, to=1, col=(4),
      add = T, lty=6, type = "l",  lwd = 2.0)
legend("topleft", c(expression(paste(plain(a), " = 0.8, ", 
                                     plain(b), " = 2.5, ",
                                     plain(alpha), " = 2.0, ",
                                     plain(lambda), " = 0.5 ")),
                    c(expression(paste(plain(a), " = 1.5, ", 
                                       plain(b), " = 2.0, ",
                                       plain(alpha), " = 2.0, ",
                                       plain(lambda), " = 0.4 ")),
                      c(expression(paste(plain(a), " = 1.4, ", 
                                         plain(b), " = 0.8, ",
                                         plain(alpha), " = 3.0, ",
                                         plain(lambda), " = 0.5 ")),
                        c(expression(paste(plain(a), " = 2.0, ", 
                                           plain(b), " = 1.5, ",
                                           plain(alpha), " = 1.8, ",
                                           plain(lambda), " = 1.5 ")))))),
       col = c(1,2,3,4), lty= c(1,2,3,6), lwd = c(3,3,6,3), 
       bty="n", cex = 2.2)
dev.off()

postscript(file = "Bp.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
a=0.8; b=2.5; alpha=2; lambda=.5
curve(Bp,from=0.00095, to=1, 
      add = FALSE, lty=1, cex.lab=2.5, type = "l", xlab = "p",
      ylab = expression("B(p)"), cex.axis=2, ylim =c(0,1.22), lwd = 2.0)
a=1.5; b=2; alpha=2; lambda=.4
curve(Bp,from=0.00095, to=1, col=(2),
        add = T, lty=2, type = "l",  lwd = 2.0)
a=1.4; b=.8; alpha=3; lambda=.5
curve(Bp,from=0.00095, to=1, col=(3),
      add = T, lty=3, type = "l",  lwd = 6.0)
a=2; b=1.5; alpha=1.8; lambda=1.5
curve(Bp,from=0.00095, to=1, col=(4),
      add = T, lty=6, type = "l",  lwd = 2.0)
legend("topleft", c(expression(paste(plain(a), " = 0.8, ", 
                                      plain(b), " = 2.5, ",
                                      plain(alpha), " = 2.0, ",
                                      plain(lambda), " = 0.5 ")),
                     c(expression(paste(plain(a), " = 1.5, ", 
                                        plain(b), " = 2.0, ",
                                        plain(alpha), " = 2.0, ",
                                        plain(lambda), " = 0.4 ")),
                       c(expression(paste(plain(a), " = 1.4, ", 
                                          plain(b), " = 0.8, ",
                                          plain(alpha), " = 3.0, ",
                                          plain(lambda), " = 0.5 ")),
                         c(expression(paste(plain(a), " = 2.0, ", 
                                            plain(b), " = 1.5, ",
                                            plain(alpha), " = 1.8, ",
                                            plain(lambda), " = 1.5 ")))))),
       col = c(1,2,3,4), lty= c(1,2,3,6), lwd = c(3,3,6,3), 
       bty="n", cex = 2.2)
dev.off()

