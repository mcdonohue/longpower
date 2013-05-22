\name{power.mmrm}
\alias{power.mmrm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear mixed model sample size calculations.
}
\description{
This function performs the sample size calculation for a mixed model of repeated measures with general correlation structure. See Lu, Luo, & Chen (2008) for parameter definitions and other details. This function executes Formula (3) on page 4.
}
\usage{
power.mmrm(N = NULL, Ra, ra, sigmaa, Rb = NULL, rb = NULL, sigmab = NULL, lambda = 1,
    delta = NULL, sig.level = 0.05, power = NULL, 
    alternative = c("two.sided", "one.sided"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
  total sample size
}
  \item{Ra}{
  correlation matrix for group a
}
  \item{ra}{
  retention in group a
}
  \item{sigmaa}{
  standard deviation of observation of interest in group a
}
  \item{Rb}{
  correlation matrix for group a
}
  \item{rb}{
  retention in group b
}
  \item{sigmab}{
  standard deviation of observation of interest in group b. If NULL, \code{sigmab} is
  assumed same as \code{sigmaa}. If not NULL, \code{sigmaa} and \code{sigmab} are averaged.
}
  \item{lambda}{
  allocation ratio
}
  \item{delta}{
  effect size
}
  \item{sig.level}{
  type one error
}
  \item{power}{
  power
}
  \item{alternative}{
  one- or two-sided test
}
}
\details{
See Lu, Luo, & Chen (2008).
}
\value{
The number of subject required per arm to attain the specified \code{power} given \code{sig.level} and the other parameter estimates.
}
\references{
Lu, K., Luo, X., Chen, P.-Y. (2008) Sample size estimation for repeated measures analysis in randomized clinical trials with missing data. \emph{International Journal of Biostatistics}, 4, (1)
}
\author{
Michael C. Donohue
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{power.mmrm.ar1}}, \code{\link{lmmpower}}, \code{\link{diggle.linear.power}}
}
\examples{
# reproduce Table 1 from Lu, Luo, & Chen (2008)
phi1 <- c(rep(1, 6), 2, 2)
phi2 <- c(1, 1, rep(2, 6))
lambda <- c(1, 2, sqrt(1/2), 1/2, 1, 2, 1, 2)
ztest <- ttest1 <- c()
for(i in 1:8){
  Na <- (phi1[i] + lambda[i] * phi2[i])*(qnorm(0.05/2) + qnorm(1-0.90))^2*(0.5^-2)
  Nb <- Na/lambda[i]
  ztest <- c(ztest, Na + Nb)
  v <- Na + Nb - 2
  Na <- (phi1[i] + lambda[i] * phi2[i])*(qt(0.05/2, df = v) + qt(1-0.90, df = v))^2*(0.5^-2)
  Nb <- Na/lambda[i]
  ttest1 <- c(ttest1, Na + Nb)
}
data.frame(phi1, phi2, lambda, ztest, ttest1)

Ra <- matrix(0.25, nrow = 4, ncol = 4)
diag(Ra) <- 1
ra <- c(1, 0.90, 0.80, 0.70)
sigmaa <- 1
power.mmrm(Ra = Ra, ra = ra, sigmaa = sigmaa, delta = 0.5, power = 0.80)
power.mmrm(N = 174, Ra = Ra, ra = ra, sigmaa = sigmaa, delta = 0.5)
power.mmrm(N = 174, Ra = Ra, ra = ra, sigmaa = sigmaa, power = 0.80)

power.mmrm(Ra = Ra, ra = ra, sigmaa = sigmaa, delta = 0.5, power = 0.80, lambda = 2)
power.mmrm(N = 174, Ra = Ra, ra = ra, sigmaa = sigmaa, delta = 0.5, lambda = 2)
power.mmrm(N = 174, Ra = Ra, ra = ra, sigmaa = sigmaa, power = 0.80, lambda = 2)
}
\keyword{ power }
\keyword{ sample size }% __ONLY ONE__ keyword per line
\keyword{ mixed effects }
\keyword{ random effects }
