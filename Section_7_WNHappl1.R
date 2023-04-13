library(AdequacyModel)
source("Section_7_functions.R")

## Proschan (1976)
x <-  c(194, 15, 41, 29, 33, 181, 413, 14, 58, 37, 100, 65, 9, 169, 447, 184, 36, 201, 118, 34, 31, 18, 
        18, 67, 57, 62, 7, 22, 34, 90, 10, 60, 186, 61, 49, 14, 24, 56, 20, 79, 84, 44, 59, 29, 118, 25,
        156, 310, 76, 26, 44, 23, 62, 130, 208, 70, 101, 208, 74, 57, 48, 29, 502, 12, 70, 21, 29, 386, 
        59, 27, 153, 26, 326, 55, 320, 56, 104, 220, 239, 47, 246, 176, 182, 33, 15, 104, 35, 23, 261, 
        87, 7, 120, 14, 62, 47, 225, 71, 246, 21, 42, 20, 5, 12, 120, 11, 3, 14, 71, 11, 14, 11, 16, 90, 
        1, 16, 52, 95, 97, 51, 11, 4, 141, 18, 142, 68, 77, 80, 1, 16, 106, 206, 82, 54, 31, 216, 46, 
        111, 39, 63, 18, 191, 18, 163, 24, 50, 44, 102, 72, 22, 39, 3, 15, 197, 188, 79, 88, 46, 5, 5, 
        36, 22, 139, 210, 97, 30, 23, 13, 14, 359, 9, 12, 270, 603, 3, 104, 2, 438, 50, 254, 5, 283, 35, 12, 
        130, 493, 487, 18, 100, 7, 98, 5, 85, 91, 43, 230, 3, 130, 102, 209, 14, 57, 54, 32, 67, 59, 134, 152, 
        27, 14, 230, 66, 61, 34)

set.seed(2000)
# Weibull Nadarajah-Haghighi
mod1 = goodness.fit(pdf=f_WNH, cdf=F_WNH, 
                    starts = c(6,3,.1,3), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)

# EGNH
mod2 = goodness.fit(pdf=f_EGNH, cdf=F_EGNH, 
                    starts = c(1,1,1,1), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)

# Kw Nadarajah-Haghighi
mod3 = goodness.fit(pdf=f_KwNH, cdf=F_KwNH, starts = c(1.7,.2,.7,.1), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)

# Beta NH
mod4 = goodness.fit(pdf=f_BNH, cdf=F_BNH, starts = c(1,1,1,1), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)

# Kw weibull
mod5 = goodness.fit(pdf=f_kww, cdf=F_kww, starts = c(1,1,1,1), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)

# Generalized power weibull
mod6 = goodness.fit(pdf=fdp_n_h, cdf=fda_n_h, starts = c(1,2,2), data = x,
                            method="S", domain=c(0,Inf),mle=NULL)
  
#exp-weibull
mod7 = goodness.fit(pdf=fdp_expweibull, cdf=fda_expweibull,
                    starts=c(.1,.2,1), data=x, method="S",
                    domain=c(0,Inf),mle=NULL)

# Exponentiated Nadarajah-Haghighi
mod8 = goodness.fit(pdf=fdp_enh, cdf=fda_enh, 
                    starts = c(1,1,1), data = x,
                    method="S", domain=c(0,Inf),mle=NULL)

# Nadarajah-Haghighi
mod9 = goodness.fit(pdf=fdp_nh, cdf=fda_nh, 
                     starts = c(1,1), data = x,
                     method="B", domain=c(0,Inf),mle=NULL)

# Weibull
mod10 = goodness.fit(pdf=fdp_we, cdf=fda_we, starts=c(1,.7), data=x,
                    method="S", domain=c(0,Inf),mle=NULL)


#Weibull exponential
set.seed(000)
mod11 = goodness.fit(pdf=f_WExp, cdf=F_WExp, starts=c(1,1,1), data=x,
                     method="S", domain=c(0,Inf),mle=NULL)
# Gompertz
set.seed(2000)
mod12 = goodness.fit(pdf=f_Go, cdf=F_Go, starts=c(2,6), data=x,
             method="S", domain=c(0,Inf),mle=NULL)

#gamma nadarajah-haghighi
set.seed(2000)
mod5 = goodness.fit(pdf=f_ZBNH, cdf=F_ZBNH, starts=c(1,.5,.1), data=x,
                     method="S", domain=c(0,Inf),mle=NULL)

estim<-matrix(NA,nrow=24,ncol = 7) # matriz que guardar?a todos os resultados 
rownames(estim)<-c("WNH", "", "EGNH", "", "KwNH", "", "BNH", "","GNH", "",
                   "GPW", "","EW","", "ENH", "","NH","","Weibull","",
                   "Wexp","","Gompertz","")
estim[1,] <- format(c(mod1$mle,mod1$KS$statistic,mod1$A,mod1$W), digits = 4)
estim[2,] <- c(format(mod1$Erro, digits = 4),rep("",3))
estim[3,] <- format(c(mod2$mle,mod2$KS$statistic,mod2$A,mod2$W), digits = 4)
estim[4,] <- c(format(mod2$Erro, digits = 4),rep("",3))
estim[5,] <- c(format(mod3$mle, digits = 4),mod3$KS$statistic,mod3$A,mod3$W)
estim[6,] <- c(format(mod3$Erro, digits = 4),rep("",3))
estim[7,] <- c(format(mod4$mle, digits = 4),mod4$KS$statistic,mod4$A,mod4$W)
estim[8,] <- c(format(mod4$Erro, digits = 4),rep("",3))
estim[9,] <- c(format(mod5$mle, digits = 4),"",mod5$KS$statistic,mod5$A,mod5$W)
estim[10,] <- c(format(mod5$Erro, digits = 4),"",rep("",3))
estim[11,] <- c(format(c(abs(mod6$mle)), digits = 4),"",mod6$KS$statistic,mod6$A,mod6$W)
estim[12,] <- c(format(c(abs(mod6$Erro)), digits = 4),"",rep("",3))
estim[13,] <- c(format(mod7$mle, digits = 4),"",mod7$KS$statistic,mod7$A,mod7$W)
estim[14,] <- c(format(mod7$Erro, digits = 4),"",rep("",3))
estim[15,] <- c(format(mod8$mle, digits = 4),"",mod8$KS$statistic,mod8$A,mod8$W)
estim[16,] <- c(format(mod8$Erro, digits = 4),"",rep("",3))
estim[17,] <- c(format(mod9$mle, digits = 4),"","",mod9$KS$statistic,mod9$A,mod9$W)
estim[18,] <- c(format(mod9$Erro, digits = 4),"","",rep("",3))
estim[19,] <- c(format(mod10$mle, digits = 4),"","",mod10$KS$statistic,mod10$A,mod10$W)
estim[20,] <- c(format(mod10$Erro, digits = 4),"","",rep("",3))
estim[21,] <- c(format(mod11$mle, digits = 4),"",mod11$KS$statistic,mod11$A,mod11$W)
estim[22,] <- c(format(mod11$Erro, digits = 4),"",rep("",3))
estim[23,] <- c(format(mod12$mle, digits = 4),"","",mod12$KS$statistic,mod12$A,mod12$W)
estim[24,] <- c(format(mod12$Erro, digits = 4),"","",rep("",3))


print(estim, quote = FALSE)

par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.1)
h11<-w1<-9 
y <- seq(from = 0,to = max(x),length.out = 10000)

setEPS()
postscript(file = "WNH_hist1.eps",width = w1, height = h11,family = "Times")
{
hist(x, xlab="", ylab="", main="",col="white", freq = F,
     ylim=c(0,.0125),xlim=c(0,600),cex.lab=2.5, cex.axis=2)
lines(y,f_WNH(c(mod1$mle),y), lty=1, cex.lab=2.5, col =1, lwd = 2.0)
lines(y,f_KwNH(c(mod3$mle),y), lty=2, type = "l", col = 2, lwd = 2.0)
lines(y, fdp_enh(c(mod8$mle),y), lty=3, type = "l", col = 3, lwd = 2.0)
legend("topright",c("WNH", "KwNH","ENH"), cex = 2.5,
       col = c(1,2,3), lty= c(1,2,3), lwd = c(3,3,3), 
       pt.bg="white", bty="n")
}
dev.off()


postscript(file = "WNH_ecdf1.eps",width = w1, height = h11,family = "Times")
{
plot(ecdf(x), main="", cex = .1, cex.lab=2.5, cex.axis=2,ylab = "")
lines(y,F_WNH(c(mod1$mle),y),lty=1,
      col=2, lwd=2,type="l")
lines(y,F_KwNH(c(mod3$mle),y),lty=2, col=3, #col="#787878",
      lwd=2,type="l")
lines(y, fda_enh(c(mod8$mle),y),lty=3,col=4,
      lwd=2,type="l")
legend("bottomright",c("WNH", "KwNH","ENH"), cex = 2.5,
       col = c(1,2,3), lty= c(1,2,3), lwd = c(3,3,3), 
       pt.bg="white", bty="n")
}
dev.off()


