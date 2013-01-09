\name{power.mmrm.ar1}
\alias{power.mmrm.ar1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear mixed model sample size calculations.
}
\description{
This function performs the sample size calculation for a mixed model of repeated measures with AR(1) correlation structure. See Lu, Luo, & Chen (2008) for parameter definitions and other details.
}
\usage{
power.mmrm.ar1(rho, ra, sigmaa, rb = NULL, sigmab = NULL, 
  lambda = 1, times = 1:length(ra),
  delta = NULL, sig.level = 0.05, power = NULL, 
  alternative = c("two.sided", "one.sided"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho}{
  AR(1) correlation parameter
}
  \item{ra}{
  retention in group a
}
  \item{sigmaa}{
  standard deviation of observation of interest in group a
}
  \item{rb}{
  retention in group a (assumed same as \code{ra} if left blank)
}
  \item{sigmab}{
  standard deviation of observation of interest in group b. If NULL, \code{sigmab} is
  assumed same as \code{sigmaa}. If not NULL, \code{sigmaa} and \code{sigmab} are averaged.
}
  \item{lambda}{
  allocation ratio
}
  \item{times}{
  observation times
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
Lu, K., Luo, X., & Chen, P.-Y. (July 14, 2008). Sample size estimation for repeated measures analysis in randomized clinical trials with missing data. International Journal of Biostatistics, 4, (1)
}
\author{
Michael C. Donohue
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{power.mmrm}}, \code{\link{lmmpower}}, \code{\link{diggle.linear.power}}
}
\examples{
# reproduce Table 2 from Lu, Luo, & Chen (2008)
tab <- c()
for(J in c(2,4))
for(aJ in (1:4)/10)
for(p1J in c(0, c(1, 3, 5, 7, 9)/10)){
  rJ <- 1-aJ
  r <- seq(1, rJ, length = J)
  # p1J = p^(J-1)
  tab <- c(tab, power.mmrm.ar1(rho = p1J^(1/(J-1)), ra = r, sigmaa = 1, 
    lambda = 1, times = 1:J,
    delta = 1, sig.level = 0.05, power = 0.80)$phi1)
}
matrix(tab, ncol = 6, byrow = TRUE)

# approximate simulation results from Table 5 from Lu, Luo, & Chen (2008)
ra <- c(100, 76, 63, 52)/100
rb <- c(100, 87, 81, 78)/100

power.mmrm.ar1(rho=0.6, ra=ra, sigmaa=1, rb = rb, lambda = sqrt(1.25/1.75), power = 0.904,
    delta = 0.9)
power.mmrm.ar1(rho=0.6, ra=ra, sigmaa=1, rb = rb, lambda = 1.25/1.75, power = 0.910,
    delta = 0.9)
power.mmrm.ar1(rho=0.6, ra=ra, sigmaa=1, rb = rb, lambda = 1, power = 0.903,
    delta = 0.9)
power.mmrm.ar1(rho=0.6, ra=ra, sigmaa=1, rb = rb, lambda = 2, power = 0.904,
    delta = 0.9)
}
\keyword{ power }
\keyword{ sample size }% __ONLY ONE__ keyword per line
\keyword{ mixed effects }
\keyword{ random effects }
