\name{edland.linear.power}
\alias{edland.linear.power}
\title{Linear mixed model sample size calculations.}
\usage{
  edland.linear.power(n = NULL, delta = NULL, t = NULL,
    sig2.s = 0, sig2.e = 1, sig.level = 0.05, power = NULL,
    alternative = c("two.sided", "one.sided"))
}
\arguments{
  \item{n}{sample size per group}

  \item{delta}{group difference in slopes}

  \item{t}{the observation times}

  \item{sig2.s}{variance of random slope}

  \item{sig2.e}{residual variance}

  \item{sig.level}{type one error}

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
  linear mixed model with random slope.
}
\details{
  This function will also provide sample size estimates for
  linear mixed models with random intercept only simply by
  setting \code{sig2.s = 0}
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
        ceiling(edland.linear.power(
          delta=0.5,
          t=t,
          sig2.e=sigma2*(1-rho),
          alternative="one.sided",
          power=0.80)$n)}))
colnames(tab) = paste("sigma2 =", sigma2)
rownames(tab) = paste("rho =", rho)
tab

# An Alzheimer's Disease example using ADAS-cog pilot estimates
t = seq(0,1.5,0.25)
n = length(t)

edland.linear.power(delta=1.5, t=t, sig2.s = 24, sig2.e = 10, sig.level=0.05, power = 0.80)
}
\author{
  Michael C. Donohue, Steven D. Edland
}
\references{
  Edland, S.D. (2009) Which MRI measure is best for
  Alzheimer's disease prevention trials: Statistical
  considerations of power and sample size. \emph{Joint Stat
  Meeting Proceedings}. 4996-4999.
}
\seealso{
  \code{\link{lmmpower}},
  \code{\link{diggle.linear.power}},
  \code{\link{liu.liang.linear.power}}
}
\keyword{effects}
\keyword{mixed}
\keyword{power}
\keyword{random}
\keyword{sample}
\keyword{size}

