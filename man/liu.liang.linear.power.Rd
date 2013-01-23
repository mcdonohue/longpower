\name{liu.liang.linear.power}
\alias{liu.liang.linear.power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear mixed model sample size calculations from Liu & Liang (1997).
}
\description{
This function performs the sample size calculation for a linear mixed model. See Liu and Liang (1997) for parameter definitions and other details.
}
\usage{
liu.liang.linear.power(N = NULL, delta = NULL, u = NULL, v = NULL, sigma2 =
     1, R = NULL, R.list = NULL, sig.level = 0.05, power =
     NULL, Pi = rep(1/length(u), length(u)), alternative =
     c("two.sided", "one.sided"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
  The total sample size. This formula can accommodate unbalanced group allocation via \code{Pi}. See Liu and Liang (1997) for more details
}
  \item{delta}{
  group difference (possibly a vector of differences)
}
  \item{u}{
%%     ~~Describe \code{u} here~~
  a list of covariate vectors or matrices associated with the parameter of interest
}
  \item{v}{
%%     ~~Describe \code{v} here~~
  a respective list of covariate vectors or matrices associated with the nuisance parameter
}
  \item{sigma2}{
  the error variance
%%     ~~Describe \code{sigma2} here~~
}
  \item{R}{
  the variance-covariance matrix for the repeated measures
%%     ~~Describe \code{R} here~~
}
  \item{R.list}{
  a list of variance-covariance matrices for the repeated measures, if assumed different in two groups
%%     ~~Describe \code{R} here~~
}
  \item{sig.level}{
  type one error
%%     ~~Describe \code{alpha} here~~
}
  \item{power}{
  power
%%     ~~Describe \code{power} here~~
}
  \item{Pi}{
  the proportion of covariates of each type
%%     ~~Describe \code{Pi} here~~
}
  \item{alternative}{
  one- or two-sided test
%%     ~~Describe \code{Pi} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The parameters \code{u}, \code{v}, and \code{Pi} are expected to be the same length and sorted with respect to each other. See Liu and Liang (1997) and package vignette for more details.
}
%\value{
%The number of subject required per arm to attain the specified \code{power} given \code{alpha} and the other parameter estimates.
%}
\references{
Liu, G., & Liang, K. Y. (1997). Sample size calculations for studies with correlated observations. \emph{Biometrics}, 53(3), 937-47.
}
%\author{
%Michael C. Donohue
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{lmmpower}}
}
\examples{
\dontrun{
browseVignettes(package = "longpower")
}
# Reproduces the table on page 29 of Diggle et al
n = 3
t = c(0,2,5)
u = list(u1 = t, u2 = rep(0,n))
v = list(v1 = cbind(1,1,rep(0,n)),
         v2 = cbind(1,0,t))         
rho = c(0.2, 0.5, 0.8)
sigma2 = c(100, 200, 300)
tab = outer(rho, sigma2, 
      Vectorize(function(rho, sigma2){
        ceiling(liu.liang.linear.power(
          delta=0.5, u=u, v=v,
          sigma2=sigma2,
          R=rho, alternative="one.sided",
          power=0.80)$N/2)}))
colnames(tab) = paste("sigma2 =", sigma2)
rownames(tab) = paste("rho =", rho)
tab

# An Alzheimer's Disease example using ADAS-cog pilot estimates
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

liu.liang.linear.power(delta=1.5, u=u, v=v, R=R, sig.level=0.05, power=0.80)
liu.liang.linear.power(N=416, u=u, v=v, R=R, sig.level=0.05, power=0.80)
liu.liang.linear.power(N=416, delta = 1.5, u=u, v=v, R=R, sig.level=0.05)
liu.liang.linear.power(N=416, delta = 1.5, u=u, v=v, R=R, power=0.80, sig.level = NULL)

# Reproduces total sample sizes, m, of Table 1 of Liu and Liang 1997
tab1 <- data.frame(cbind(
  n = c(rep(4, 4), rep(2, 4), 1),
  rho = c(0.0, 0.3, 0.5, 0.8)))
u = list(u1 = 1, u2 = 1) # intercept
v = list(v1 = 1, # treatment
       v2 = 0) # control       
m <- c()
for(i in 1:nrow(tab1)){
  R <- matrix(tab1$rho[i], nrow = tab1$n[i], ncol = tab1$n[i])
  diag(R) <- 1
  m <- c(m, liu.liang.linear.power(
    delta=0.5,
    u = list(u1 = rep(1, tab1$n[i]), u2 = rep(1, tab1$n[i])), # intercept
    v = list(v1 = rep(1, tab1$n[i]), # treatment
             v2 = rep(0, tab1$n[i])), # control       
    sigma2=1,
    R=R, alternative="two.sided",
    power=0.90)$N)
}
cbind(tab1, m)

# Reproduces total sample sizes, m, of Table 3.a. of Liu and Liang 1997
# with unbalanced design
tab3 <- data.frame(cbind(
  rho = rep(c(0.0, 0.3, 0.5, 0.8), 2),
  pi1 = c(rep(0.8, 4), rep(0.2, 4))))
u = list(u1 = 1, u2 = 1) # intercept
v = list(v1 = 1, # treatment
       v2 = 0) # control
m <- c()
for(i in 1:nrow(tab3)){
  R <- matrix(tab3$rho[i], nrow = 4, ncol = 4)
  diag(R) <- 1
  m <- c(m, ceiling(liu.liang.linear.power(
    delta=0.5,
    u = list(u1 = rep(1, 4), u2 = rep(1, 4)), # intercept
    v = list(v1 = rep(1, 4), # treatment
             v2 = rep(0, 4)), # control       
    sigma2=1,
    Pi = c(tab3$pi1[i], 1-tab3$pi1[i]),
    R=R, alternative="two.sided",
    power=0.90)$N))
}
cbind(tab3, m)



}
\keyword{ power }
\keyword{ sample size }% __ONLY ONE__ keyword per line
\keyword{ mixed effects }
\keyword{ random effects }
