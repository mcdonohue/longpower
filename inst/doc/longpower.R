### R code from vignette source 'longpower.Rnw'

###################################################
### code chunk number 1: longpower.Rnw:75-77
###################################################
options(prompt=" ", continue="  ")
#library(mvtnorm)


###################################################
### code chunk number 2: longpower.Rnw:79-80
###################################################
library(longpower)


###################################################
### code chunk number 3: longpower.Rnw:82-83
###################################################
liu.liang.linear.power


###################################################
### code chunk number 4: longpower.Rnw:89-108
###################################################
n = 3
t = c(0,2,5)
u = list(u1 = t, u2 = rep(0,n))
v = list(v1 = cbind(1,1,rep(0,n)),
         v2 = cbind(1,0,t))         
rho = c(0.2, 0.5, 0.8)
sigma2 = c(100, 200, 300)
tab = outer(rho, sigma2, 
      Vectorize(function(rho, sigma2){
        round(diggle.linear.power(
          d=0.5,
          t=t,
          sigma2=sigma2,
          R=rho,
          alternative="one.sided",
          power=0.80)$n)}))
colnames(tab) = paste("sigma2 =", sigma2)
rownames(tab) = paste("rho =", rho)
tab


###################################################
### code chunk number 5: longpower.Rnw:112-134
###################################################
# var of random intercept
sig2.i = 55
# var of random slope
sig2.s = 24
# residual var
sig2.e = 10
# covariance of slope and intercep
cov.s.i <- 0.8*sqrt(sig2.i)*sqrt(sig2.s)

cov.t <- function(t1, t2, sig2.i, sig2.s, cov.s.i){
        sig2.i + t1*t2*sig2.s + (t1+t2)*cov.s.i 
}

t = seq(0,1.5,0.25)
n = length(t)
R = outer(t, t, function(x,y){cov.t(x,y, sig2.i, sig2.s, cov.s.i)})
R = R + diag(sig2.e, n, n)
u = list(u1 = t, u2 = rep(0,n))
v = list(v1 = cbind(1,1,rep(0,n)),
         v2 = cbind(1,0,t))         

liu.liang.linear.power(d=1.5, u=u, v=v, R=R, sig.level=0.05, alternative="two.sided", power=0.80)


###################################################
### code chunk number 6: longpower.Rnw:158-160
###################################################
x = (rbind(1,t)%*%solve(R)%*%cbind(1,t))[2,2]
x*2*(qnorm(1-0.05/2) + qnorm(0.80))^2/1.5^2


###################################################
### code chunk number 7: longpower.Rnw:182-184
###################################################
x = solve(rbind(1,t)%*%solve(R)%*%cbind(1,t))[2,2]
x*2*(qnorm(1-0.05/2) + qnorm(0.80))^2/1.5^2


###################################################
### code chunk number 8: longpower.Rnw:226-232
###################################################
X = t(v[[1]])%*%solve(R)%*%v[[1]] + 
    t(v[[2]])%*%solve(R)%*%v[[2]]

Sigma1 = ((t(u[[1]])%*%solve(R)%*%t-t(u[[1]])%*%solve(R)%*%v[[1]]%*%solve(X)%*%t(v[[1]])%*%solve(R)%*%t)/2)

(qnorm(1-0.05/2) + qnorm(0.80))^2/(Sigma1*(1.5)^2)


