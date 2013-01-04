\name{lmmpower-methods}
\docType{methods}
\alias{lmmpower-methods}
\alias{lmmpower}
\alias{lmmpower.default}
\alias{lmmpower.lme}
\alias{lmmpower.gee}
\alias{lmmpower,ANY-method}
\alias{lmmpower,mer-method}
\title{Sample size calculations for linear mixed models of rate of change based on lmer, lme, or gee "placebo" pilot estimates.
}
\description{These functions compute
        sample size for linear mixed models based on the 
        formula due to Diggle (2002) or Liu and Liang (1997). 
        These formulae are expressed
        in terms of marginal model or Generalized Estimating Equations (GEE) parameters.
        These functions translate pilot mixed 
        effect model parameters (e.g. random intercept and/or slope, fixed effects, etc.)
        into marginal model parameters so that either formula can be applied to equivalent
        affect. Pilot estimates are assumed to be from an appropriate "placebo" group and the
        parameter of interest is assumed to be the rate of change over time
        of the outcome.
}
%\section{Methods}{
%\describe{
%
%\item{\code{signature(object = "ANY")}}{ ~~describe this method here }
%
%\item{\code{signature(object = "mer")}}{ ~~describe this method here }
%}}

\usage{

lmmpower.default(object=NULL,
   parameter = 2,
   pct.change = NULL,
   delta = NULL,
   t = NULL,
   sig.level = 0.05,
   power = 0.80, 
   alternative = c("two.sided", "one.sided"),
   beta=NULL,
   beta.lim=NULL,
   delta.lim=NULL,
   sig2.i=NULL,
   sig2.s=NULL,
   sig2.e=NULL,
   cov.s.i=NULL,
   R=NULL,
   method=c("diggle", "liuliang"),
   ...)

lmmpower.lme(object,
   parameter = 2,
   pct.change = NULL,
   delta = NULL,
   t = NULL,
   sig.level = 0.05,
   power = 0.80, 
   alternative = c("two.sided", "one.sided"),
   beta=NULL,
   beta.lim=NULL,
   delta.lim=NULL,
   sig2.i=NULL,
   sig2.s=NULL,
   sig2.e=NULL,
   cov.s.i=NULL,
   ...)
   
lmmpower.gee (object,
   parameter = 2,
   pct.change = NULL,
   delta = NULL,
   t = NULL,
   sig.level = 0.05,
   power = 0.80, 
   alternative = c("two.sided", "one.sided"),
   beta=NULL,
   beta.lim=NULL,
   delta.lim=NULL,
   ...)
}

\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
an object returned by lme4 of a mixed-effects model object to placebo data assumed to have either a random intercept, or a random intercept and random effect for time (slope); and fixed effect representing the rate of change in a placebo group.
}
  \item{parameter}{
%%     ~~Describe \code{parameter} here~~
the name or position of the rate of change parameter of interest, e.g. (\code{"time"}, \code{"t"}, or \code{2} if it is the second specified fixed effect).
}
  \item{pct.change}{
%%     ~~Describe \code{pct.change} here~~
the percent change in the pilot estimate of the parameter of interest (\code{beta}, the placebo/null effect) 
}
  \item{delta}{
%%     ~~Describe \code{d} here~~
the change in the pilot estimate of the parameter of interest, computed from \code{pct.change} if left missing.
}
  \item{t}{
%%     ~~Describe \code{t} here~~
vector of time points
}
  \item{sig.level}{
%%     ~~Describe \code{sig.level} here~~
Type I error
}
  \item{power}{
%%     ~~Describe \code{power} here~~
power
}
  \item{alternative}{
%%     ~~Describe \code{alternative} here~~
\code{"two.sided"} or \code{"one.sided"}
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
pilot estimate of the placebo effect (slope or rate of change in the outcome)
}
  \item{beta.lim}{
%%     ~~Describe \code{beta.lim} here~~
95\% confidence limits of the pilot estimate of beta
}
  \item{delta.lim}{
%%     ~~Describe \code{d.lim} here~~
95\% confidence limits of the effect size
}
  \item{sig2.i}{
%%     ~~Describe \code{sig2.i} here~~
pilot estimate of variance of random intercept
}
  \item{sig2.s}{
%%     ~~Describe \code{sig2.s} here~~
pilot estimate of variance of random slope
}
  \item{sig2.e}{
%%     ~~Describe \code{sig2.e} here~~
pilot estimate of residual variance
}
  \item{cov.s.i}{
%%     ~~Describe \code{cov.s.i} here~~
pilot estimate of covariance of random slope and intercept
}
  \item{R}{
%%     ~~Describe \code{cov.s.i} here~~
pilot estimate of a marginal model working correlation matrix
}
  \item{method}{
%%     ~~Describe \code{\dots} here~~
the formula to use. Defaults to \code{"diggle"} for Diggle et al (2002). Alternatively
\code{"liuliang"} can be selected for Liu & Liang (1997).
}
  \item{...}{
%%     ~~Describe \code{\dots} here~~
other arguments}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
Any parameters not explicitly stated are extracted from the fitted \code{object}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of class \code{power.htest} giving the pasm
}
\references{
Diggle PJ, Heagerty PJ, Liang K-, Zeger SL. Analysis of longitudinal data. Second Edition. Oxford Statistical Science Serires. 2002.

Liu, G., & Liang, K. Y. (1997). Sample size calculations for studies with correlated observations. Biometrics, 53(3), 937-47.
}
\author{
Michael C. Donohue
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{liu.liang.linear.power}}
\code{\link{diggle.linear.power}}
}
\examples{
lmmpower(delta=1.5, t = seq(0,1.5,0.25),
	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24))

lmmpower(beta = 5, pct.change = 0.30, t = seq(0,1.5,0.25),
	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24))

\dontrun{
library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
lmmpower(fm1, pct.change = 0.30, t = seq(0,9,1))

library(nlme)
fm2 <- lme(Reaction ~ Days, random=~Days|Subject, sleepstudy)
lmmpower(fm2, pct.change = 0.30, t = seq(0,9,1))

# random intercept only
fm3 <- lme(Reaction ~ Days, random=~1|Subject, sleepstudy)
lmmpower(fm3, pct.change = 0.30, t = seq(0,9,1))

library(gee)
fm4 <- gee(Reaction ~ Days, id = Subject,
            data = sleepstudy,
            corstr = "exchangeable")
lmmpower(fm4, pct.change = 0.30, t = seq(0,9,1))

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ power }
\keyword{ sample size }% __ONLY ONE__ keyword per line
\keyword{ mixed effects }
\keyword{ random effects }
\keyword{ marginal model }
\keyword{methods}
