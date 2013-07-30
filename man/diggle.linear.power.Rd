\name{diggle.linear.power}
\alias{diggle.linear.power}
\title{Linear mixed model sample size calculations.}
\usage{
  diggle.linear.power(n = NULL, delta = NULL, t = NULL,
    sigma2 = 1, R = NULL, sig.level = 0.05, power = NULL,
    alternative = c("two.sided", "one.sided"))
}
\arguments{
  \item{n}{sample size per group}

  \item{delta}{group difference in slopes}

  \item{t}{the observation times}

  \item{sigma2}{the marginal model (GEE) scale parameter}

  \item{R}{the working correlation matrix (or
  variance-covariance matrix if \code{sigma2} is 1). If
  \code{R} is a scalar, an exchangeable working correlation
  matrix will be assumed.}

  \item{sig.level}{Type I error}

  \item{power}{power}

  \item{alternative}{one- or two-sided test}
}
\value{
  The number of subject required per arm to attain the
  specified \code{power} given \code{sig.level} and the
  other parameter estimates.
}
\description{
  This function performs the sample size calculation for a
  linear mixed model. See Diggle et al (2002) for parameter
  definitions and other details.
}
\details{
  The parameters \code{u}, \code{v}, and \code{Pi} are
  expected to be the same length and sorted with respect to
  each other. See Diggle, et al (1997) and package vignette
  for more details.
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
        ceiling(diggle.linear.power(
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
\author{
  Michael C. Donohue, Steven D. Edland
}
\references{
  Diggle P.J., Heagerty P.J., Liang K., Zeger S.L. (2002)
  \emph{Analysis of longitudinal data}. Second Edition.
  Oxford Statistical Science Series.
}
\seealso{
  \code{\link{lmmpower}}, \code{\link{diggle.linear.power}}
}
\keyword{effects}
\keyword{mixed}
\keyword{power}
\keyword{random}
\keyword{sample}
\keyword{size}

