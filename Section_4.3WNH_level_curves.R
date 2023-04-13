rm(list = ls())
setwd("~/Insync/fernando.p.ramirez@ufsm.br/Google Drive - Shared with me/UFSM/Papers 2017/Weibull NH/RCE/figure")

Q <- function(u,a,b,alpha,lambda) {
    (1/lambda)*((1+log(1+((-log(1-u)/a)^(1/b))))^(1/alpha)-1)
}
B <-function(a,b){
    return((Q(3/4,a,b,alpha,lambda)-2*Q(1/2,a,b,alpha,lambda)+Q(1/4,a,b,alpha,lambda))/(Q(3/4,a,b,alpha,lambda)-Q(1/4,a,b,alpha,lambda)))
}
M <-function(a,b){
  return((Q(7/8,a,b,alpha,lambda)-Q(5/8,a,b,alpha,lambda)+Q(3/8,a,b,alpha,lambda)-Q(1/8,a,b,alpha,lambda))/(Q(3/4,a,b,alpha,lambda)-Q(1/4,a,b,alpha,lambda)))
}
####################
# Plots of B and M #
####################
h11 <- w1 <- 9 # width and height of the plots
setEPS()
alpha=2
lambda=.5

# a)
postscript(file = "B_skew_a(CL).eps",width = w1, height = h11,family = "Times")
a_range <- seq(0.02, 2.5, length.out = 100)
b_range <- seq(0.01, 1, length.out = 100)
B_vals <- outer(a_range, b_range, B)
contour(a_range, b_range, B_vals, 
       xlab = "a", ylab = "b",  cex.lab = 2.5, cex.axis = 2, labcex = 1.5)
dev.off()

#c)
postscript(file = "M_kurt_a(CL).eps",width = w1, height = h11,family = "Times")
b_range <- seq(0.15, 0.6, length.out = 1000)
a_range <- seq(0.25, 3, length.out = 1000)
M_vals <- outer(a_range, b_range, M)
contour(a_range, b_range, M_vals, 
        xlab = "a", ylab = "b",  cex.lab = 2.5, cex.axis = 2, labcex = 1.5)
dev.off()

#b)
alpha=0.8
lambda=1.5
postscript(file = "B_skew_b(CL).eps",width = w1, height = h11,family = "Times")
a_range <- seq(0.15, 1.5, length.out = 100)
b_range <- seq(0.01, 2, length.out = 100)
B_vals <- outer(a_range, b_range, B)
contour(a_range, b_range, B_vals, 
        xlab = "a", ylab = "b",  cex.lab = 2.5, cex.axis = 2, labcex = 1.5)
dev.off()

# d)
postscript(file = "M_kurt_b(CL).eps",width = w1, height = h11,family = "Times")
b_range <- seq(0.15, 0.5, length.out = 100)
a_range <- seq(0.25, 4, length.out = 100)
M_vals <- outer(a_range, b_range, M)
contour(a_range, b_range, M_vals, 
        xlab = "a", ylab = "b",  cex.lab = 2.5, cex.axis = 2, labcex = 1.5)
dev.off()

