####################################
##   Weibull Nadarajah-Haghighi   ##
####################################
## WNH cdf
F_WNH<- function(x){
  1-exp(-a*(exp((1+lambda*x)^alpha-1)-1)^b)
}

## WNH pdf
f_WNH <- function(x){
  a*b*lambda*alpha*(1+lambda*x)^(alpha-1)*(1-exp(1-(1+lambda*x)^alpha))^(b-1)*
    exp(-b*(1-(1+lambda*x)^alpha)-a*(exp((1+lambda*x)^alpha-1)-1)^b)
}

# WNH hrf
h_WNH<-function(x){
  a*b*lambda*alpha*(1+lambda*x)^(alpha-1)*(1-exp(1-(1+lambda*x)^alpha))^(b-1)*
    exp(-b*(1-(1+lambda*x)^alpha))
}

# Checking the functions
a<-0.8;b<-2.5;alpha<-2;lambda<-0.5
integrate(f_WNH, 0, Inf)
x<-0.4
F_WNH(x); integrate(f_WNH, 0, x)$value
f_WNH(x)/(1-F_WNH(x)); h_WNH(x)

h11 <- w1 <- 9 # width and height of the plots
#############
# pdf plots #
#############
setEPS()
# pdf 4
postscript(file = "dens_WNH_a.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
fromx = 0.00095; tox = 1.25
yliminf = 0; ylimsup = 3
b = 2.5; alpha = 2; lambda = .5
a = .8
curve(f_WNH,from=fromx, to=tox, add = FALSE, lty=1, cex.lab=2.5, 
      type = "l", xlab = "x", ylab = expression("f(x)"), cex.axis=2,
      ylim =c(yliminf,ylimsup), col =1, lwd = 2.0)
a = 1.5
curve(f_WNH, add = TRUE, lty=2, type = "l", col = 2, lwd = 2.0)
a = 2
curve(f_WNH, add = TRUE, lty=3, type = "l", col = 3, lwd = 2.0)
a = 3
curve(f_WNH, add = TRUE, lty=6, type = "l", col = 4, lwd = 2.0)

legend("topright", c(expression(paste(plain(a), " = 0.8 ")),
                     c(expression(paste(plain(a), " = 1.5")),
                       c(expression(paste(plain(a), " = 2.0")),
                         c(expression(paste(plain(a), " = 3.0")))))),
       col = c(1,2,3,4), lty= c(1,2,3,6), lwd = c(3,3,3,3), 
       bty="n", cex = 2.5)
dev.off()
# pdf 2
postscript(file = "dens_WNH_b.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
fromx = 0.00095; tox = 1.5
yliminf = 0; ylimsup = 2.5
a = 1.5; alpha = 2; lambda = 0.4; 
b = 1.5
curve(f_WNH,from=fromx, to=tox, add = FALSE, lty=1, cex.lab=2.5, 
      type = "l", xlab = "x", ylab = expression("f(x)"), cex.axis=2,
      ylim =c(yliminf,ylimsup), col =1, lwd = 2.0)
b = 2
curve(f_WNH, add = TRUE, lty=2, col = 2, lwd = 2.0)
b = 2.5
curve(f_WNH, add = TRUE, lty=3, col = 3, lwd = 2.0)
b = 3
curve(f_WNH, add = TRUE, lty=6, col = 4, lwd = 2.0)
legend("topright", c(expression(paste(plain(b), " = 1.5 ")),
                     c(expression(paste(plain(b), " = 2.0")),
                       c(expression(paste(plain(b), " = 2.5")),
                         c(expression(paste(plain(b), " = 3.0")))))),
       col = c(1,2,3,4), lty= c(1,2,3,6), lwd = c(3,3,3,3), 
       bty="n", cex = 2.5)
dev.off()
# pdf 3
postscript(file = "dens_WNH_c.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
fromx = 0.0000095; tox = 1
yliminf = 0; ylimsup = 3.5
a = 1.4; b = .8; lambda = .5
alpha = .8
curve(f_WNH,from=fromx, to=tox, add = FALSE, lty=1, cex.lab=2.5, 
      type = "l", xlab = "x", ylab = expression("f(x)"), cex.axis=2,
      ylim =c(yliminf,ylimsup), col =1, lwd = 2.0)
alpha = 1.3
curve(f_WNH, add = TRUE, lty=2, type = "l", col = 2, lwd = 2.0)
alpha = 2
curve(f_WNH, add = TRUE, lty=3, type = "l", col = 3, lwd = 2.0)
alpha = 3
curve(f_WNH, add = TRUE, lty=6, col = 4, lwd = 2.0)
legend("topright", c(expression(paste(plain(alpha), " = 0.8 ")),
                     c(expression(paste(plain(alpha), " = 1.3")),
                       c(expression(paste(plain(alpha), " = 2.0")),
                         c(expression(paste(plain(alpha), " = 3.0")))))),
       col = c(1,2,3,4), lty= c(1,2,3,6), lwd = c(3,3,3,3), 
       bty="n", cex = 2.5)
dev.off()
# pdf 4
postscript(file = "dens_WNH_d.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
fromx = 0.00095; tox = 1
yliminf = 0; ylimsup = 5
a = 2; b = 1.5; alpha = 1.8
lambda = .5
curve(f_WNH,from=fromx, to=tox, add = FALSE, lty=1, cex.lab=2.5, 
      type = "l", xlab = "x", ylab = expression("f(x)"), cex.axis=2,
      ylim =c(yliminf,ylimsup), col =1, lwd = 2.0)
lambda = .8
curve(f_WNH, add = TRUE, lty=2, type = "l", col = 2, lwd = 2.0)
lambda = 1.2
curve(f_WNH, add = TRUE, lty=3, type = "l", col = 3, lwd = 2.0)
lambda = 1.5
curve(f_WNH, add = TRUE, lty=6, type = "l", col = 4, lwd = 2.0)
legend("topright", c(expression(paste(plain(lambda), " = 0.5 ")),
                     c(expression(paste(plain(lambda), " = 0.8")),
                       c(expression(paste(plain(lambda), " = 1.2")),
                         c(expression(paste(plain(lambda), " = 1.5")))))),
       col = c(1,2,3,4), lty= c(1,2,3,6), lwd = c(3,3,3,3), 
       bty="n", cex = 2.5)
dev.off()
#############
# hrf plots #
#############
# hrf 1
postscript(file = "WNH_hrf1.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
fromx = 0.00000095; tox = .6
yliminf = 0; ylimsup = 3
a = 1.2; b = 0.3; alpha = 0.8; lambda = .2
curve(h_WNH,from=fromx, to=tox, add = FALSE, lty=1, cex.lab=2.5, 
      type = "l", xlab = "x", ylab = expression("h(x)"), cex.axis=2,
      ylim =c(yliminf,ylimsup), col =1, lwd = 2.0)
legend("topright",  c(expression(paste("a", " = 0.5 ")),
                      c(expression(paste(plain(b), " = 0.3")),
                        c(expression(paste(plain(alpha), " = 0.8")),
                          c(expression(paste(plain(lambda), " = 0.2")))))),
       ncol = 2,lwd = c(3), bty="n", cex = 2.5)
dev.off()
# hrf 2
postscript(file = "WNH_hrf2.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
fromx = 0.00000095; tox = 10
yliminf = 0; ylimsup = 4
a = 2.5; b = 3; alpha = 0.3; lambda = 1
curve(h_WNH,from=fromx, to=tox, add = FALSE, lty=1, cex.lab=2.5, 
      type = "l", xlab = "x", ylab = expression("h(x)"), cex.axis=2,
      ylim =c(yliminf,ylimsup), col =1, lwd = 2.0)
legend("topleft",  c(expression(paste("a", " = 2.5 ")),
                      c(expression(paste(plain(b), " = 3.0")),
                        c(expression(paste(plain(alpha), " = 0.3")),
                          c(expression(paste(plain(lambda), " = 1.0")))))),
       ncol = 2,lwd = c(3), bty="n", cex = 2.5)
dev.off()
# hrf 2
postscript(file = "WNH_hrf3.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
fromx = 0.00000095; tox = 13
yliminf = 0.4; ylimsup = 0.62
a = 1.5; b = 0.8; alpha = 0.5; lambda = 0.6
curve(h_WNH,from=fromx, to=tox, add = FALSE, lty=1, cex.lab=2.5, 
      type = "l", xlab = "x", ylab = expression("h(x)"), cex.axis=2,
      ylim =c(yliminf,ylimsup), col =1, lwd = 2.0)
legend("top",  c(expression(paste("a", " = 1.5 ")),
                     c(expression(paste(plain(b), " = 0.8")),
                       c(expression(paste(plain(alpha), " = 0.5")),
                         c(expression(paste(plain(lambda), " = 0.6")))))),
       ncol = 2,lwd = c(3), bty="n", cex = 2.5)
dev.off()
# hrf 3
postscript(file = "WNH_hrf4.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
fromx = 0.00000095; tox = 1.5
yliminf = 0.8; ylimsup = 0.95
a = 2; b = 1.1; alpha = 0.3; lambda = 2
curve(h_WNH,from=fromx, to=tox, add = FALSE, lty=1, cex.lab=2.5, 
      type = "l", xlab = "x", ylab = expression("h(x)"), cex.axis=2,
      ylim =c(yliminf,ylimsup), col =1, lwd = 2.0)
legend("topright",  c(expression(paste("a", " = 2 ")),
                 c(expression(paste(plain(b), " = 1.1")),
                   c(expression(paste(plain(alpha), " = 0.3")),
                     c(expression(paste(plain(lambda), " = 2.0")))))),
       ncol = 2, lwd = c(3), bty="n", cex = 2.5)
dev.off()