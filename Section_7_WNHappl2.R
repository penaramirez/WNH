library(AdequacyModel)
source("Section_7_functions.R")

## 
x <-  c(0.046, 1.436, 2.592, 0.140, 1.492, 2.600, 0.150, 1.580, 2.670, 0.248, 1.719, 2.717,
        0.280, 1.794, 2.819, 0.313, 1.915, 2.820, 0.389, 1.920, 2.878, 0.487, 1.963, 2.950,
        0.622, 1.978, 3.003, 0.900, 2.053, 3.102, 0.952, 2.065, 3.304, 0.996, 2.117, 3.483,
        1.003, 2.137, 3.500, 1.010, 2.141, 3.622, 1.085, 2.163, 3.665, 1.092, 2.183, 3.695,
        1.152, 2.240, 4.015, 1.183, 2.341, 4.628, 1.244, 2.435, 4.806, 1.249, 2.464, 4.881,
        1.262, 2.543, 5.140)

set.seed(2000)
# Weibull Nadarajah-Haghighi
mod1 = goodness.fit(pdf=f_WNH, cdf=F_WNH, 
                    starts = c(.1,1,.8,1), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)

# EGNH
set.seed(2000)
mod2 = goodness.fit(pdf=f_EGNH, cdf=F_EGNH, 
                    starts = c(.1,7,1,.5), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)

# Kw Nadarajah-Haghighi
set.seed(2000)
mod3 = goodness.fit(pdf=f_KwNH, cdf=F_KwNH, starts = c(5,.8,3,.2), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)

# # Beta NH
# mod4 = goodness.fit(pdf=f_BNH, cdf=F_BNH, starts = c(1,1,1,1), data = x,
#                     method="S", domain=c(0,Inf),mle=NULL)
# Kw weibull
#mod4 = goodness.fit(pdf=f_kww, cdf=F_kww, starts = c(.4,.1,1,.1), data = x,
#                    method="S", domain=c(0,Inf),mle=NULL)
# Generalized power weibull
set.seed(2000)
mod4 = goodness.fit(pdf=fdp_n_h, cdf=fda_n_h, starts = c(1.2,3,1), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)

#exp-weibull
set.seed(2000)
mod5 = goodness.fit(pdf=fdp_expweibull, cdf=fda_expweibull,
                    starts=c(1,1,.4), data=x, method="S",
                    domain=c(0,Inf),mle=NULL)

# Exponentiated Nadarajah-Haghighi
set.seed(2000)
mod6 = goodness.fit(pdf=fdp_enh, cdf=fda_enh, 
                    starts = c(3,.1,.9), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)

# Nadarajah-Haghighi
set.seed(2000)
mod7 = goodness.fit(pdf=fdp_nh, cdf=fda_nh, 
                    starts = c(2,.1), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)
#weibull
set.seed(2000)
mod8 = goodness.fit(pdf=fdp_we, cdf=fda_we, starts=c(1,1), data=x,
                    method="S", domain=c(0,Inf),mle=NULL)
#weibull exponential
set.seed(2000)
mod9 = goodness.fit(pdf=f_WExp, cdf=F_WExp, starts=c(4,3,.5), data=x,
                    method="S", domain=c(0,Inf),mle=NULL)

#gompertz
set.seed(2000)
mod10 = goodness.fit(pdf=f_Go, cdf=F_Go, starts=c(.1,.5), data=x,
                     method="S", domain=c(0,Inf),mle=NULL)

#gamma nadarajah-haghighi
set.seed(2000)
mod11 = goodness.fit(pdf=f_ZBNH, cdf=F_ZBNH, starts=c(1,1,1), data=x,
                     method="S", domain=c(0,Inf),mle=NULL)


estim<-matrix(NA,nrow=22,ncol = 7) # matriz que guardar?a todos os resultados 
rownames(estim)<-c("WNH", "", "EGNH", "", "KwNH", "","GPW", "","EW","",
                   "ENH", "","NH","","Weibull","", "WExp", "","Go", "","GNH", "")
estim[1,] <- format(c(mod1$mle,mod1$KS$statistic,mod1$A,mod1$W), digits = 4)
estim[2,] <- c(format(mod1$Erro, digits = 4),rep("",3))
estim[3,] <- format(c(mod2$mle,mod2$KS$statistic,mod2$A,mod2$W), digits = 4)
estim[4,] <- c(format(mod2$Erro, digits = 4),rep("",3))
estim[5,] <- c(format(mod3$mle, digits = 4),mod3$KS$statistic,mod3$A,mod3$W)
estim[6,] <- c(format(mod3$Erro, digits = 4),rep("",3))
estim[7,] <- c(format(mod4$mle, digits = 4),"",mod4$KS$statistic,mod4$A,mod4$W)
estim[8,] <- c(format(mod4$Erro, digits = 4),"",rep("",3))
estim[9,] <- c(format(mod5$mle, digits = 4),"",mod5$KS$statistic,mod5$A,mod5$W)
estim[10,] <- c(format(mod5$Erro, digits = 4),"",rep("",3))
estim[11,] <- c(format(c(abs(mod6$mle)), digits = 4),"",mod6$KS$statistic,mod6$A,mod6$W)
estim[12,] <- c(format(c(abs(mod6$Erro)), digits = 4),"",rep("",3))
estim[13,] <- c(format(mod7$mle, digits = 4),"","",mod7$KS$statistic,mod7$A,mod7$W)
estim[14,] <- c(format(mod7$Erro, digits = 4),"","",rep("",3))
estim[15,] <- c(format(mod8$mle, digits = 4),"","",mod8$KS$statistic,mod8$A,mod8$W)
estim[16,] <- c(format(mod8$Erro, digits = 4),"","",rep("",3))
estim[17,] <- c(format(mod9$mle, digits = 4),"",mod9$KS$statistic,mod9$A,mod9$W)
estim[18,] <- c(format(mod9$Erro, digits = 4),"",rep("",3))
estim[19,] <- c(format(mod10$mle, digits = 4),"","",mod10$KS$statistic,mod10$A,mod10$W)
estim[20,] <- c(format(mod10$Erro, digits = 4),"","",rep("",3))
estim[21,] <- c(format(mod11$mle, digits = 4),"",mod11$KS$statistic,mod11$A,mod11$W)
estim[22,] <- c(format(mod11$Erro, digits = 4),"",rep("",3))


print(estim, quote = FALSE)

par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.1)
h11<-w1<-9 
y <- seq(from = 0,to = max(x),length.out = 10000)

setEPS()
postscript(file = "WNH_hist2.eps",width = w1, height = h11,family = "Times")
{
hist(x, xlab="", ylab="", main="",col="white", freq = F,
     ylim=c(0,.35),xlim=c(0,7), cex.lab=2.5, cex.axis=2)
lines(y,f_WNH(c(mod1$mle),y), lty=1, col =1, lwd = 2.0)
lines(y,fdp_expweibull(c(mod5$mle),y), lty=2, type = "l", col = 2, lwd = 2.0)
lines(y, f_Go(c(mod10$mle),y), lty=3, type = "l", col = 3, lwd = 2.0)
legend("topright",c("WNH", "EW","Gompertz"), cex = 2.5,
       col = c(1,2,3), lty= c(1,2,3), lwd = c(3,3,3), 
       pt.bg="white", bty="n")
}
dev.off()


postscript(file = "WNH_ecdf2.eps",width = w1, height = h11,family = "Times")
{
plot(ecdf(x), main="",cex = .1, cex.lab=2.5, cex.axis=2,ylab = "")
lines(y,F_WNH(c(mod1$mle),y),lty=1,
      col=2, lwd=2,type="l")
lines(y,fda_expweibull(c(mod5$mle),y),lty=2, col=3, 
      lwd=2,type="l")
lines(y, F_Go(c(mod10$mle),y),lty=3,col=4,
      lwd=2,type="l")
legend("bottomright",c("WNH", "EW","Gompertz"), cex = 2.5,
       col = c(1,2,3), lty= c(1,2,3), lwd = c(3,3,3), 
       pt.bg="white", bty="n")
}
dev.off()


