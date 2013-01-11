\name{diggle.linear.power}
\alias{diggle.linear.power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear mixed model sample size calculations.
}
\description{
This function performs the sample size calculation for a linear mixed model. See Diggle et al (2002) for parameter definitions and other details.
}
\usage{
diggle.linear.power(n=NULL, delta=NULL, t=NULL, sigma2 = 1, R=NULL, sig.level = 0.05, power = NULL, alternative = c("two.sided", "one.sided"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  sample size per group
}
  \item{delta}{
  group difference in slopes
}
  \item{t}{
  the observation times
%%     ~~Describe \code{sigma2} here~~
}
  \item{sigma2}{
  the marginal model (GEE) scale parameter
%%     ~~Describe \code{sigma2} here~~
}
  \item{R}{
  the working correlation matrix (or variance-covariance matrix if \code{sigma2}
  is 1). If \code{R} is a scalar, an exchangeable working correlation matrix
  will be assumed.
%%     ~~Describe \code{R} here~~
}
  \item{sig.level}{
  type one error
%%     ~~Describe \code{sig.level} here~~
}
  \item{power}{
  power
%%     ~~Describe \code{power} here~~
}
  \item{alternative}{
  one- or two-sided test
%%     ~~Describe \code{Pi} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The parameters \code{u}, \code{v}, and \code{Pi} are expected to be the same length and sorted with respect to each other. See Diggle, et al (1997) and package vignette for more details.
}
\value{
The number of subject required per arm to attain the specified \code{power} given \code{sig.level} and the other parameter estimates.
}
\references{
Diggle PJ, Heagerty PJ, Liang K, Zeger SL. Analysis of longitudinal data. Second Edition. Oxford Statistical Science Serires. 2002.
}
\author{
Michael C. Donohue, Steven D. Edland
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{lmmpower}}, \code{\link{diggle.linear.power}}
}
\examples{
\dontrun{
browseVignettes(package = "longpower")
}

# Reproduces the table on page 29 of Diggle et al
n = 3
t = c(0,2,5)
rho = c(0.2, 0.5, 0.8)
sigma2 = c(100, 200, 300)
tab = outer(rho, sigma2, 
      Vectorize(function(rho, sigma2){
        round(diggle.linear.power(
          delta=0.5,
          t=t,
          sigma2=sigma2,
          R=rho,
          alternative="one.sided",
          power = 0.80)$n)}))
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

diggle.linear.power(d=1.5, t=t, R=R, sig.level=0.05, power=0.80)
}
\keyword{ power }
\keyword{ sample size }% __ONLY ONE__ keyword per line
\keyword{ mixed effects }
\keyword{ random effects }
