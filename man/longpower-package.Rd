\name{longpower-package}
\alias{longpower-package}
\docType{package}
\title{Sample size calculations for longitudinal data}
\description{The longpower package contains functions for computing
        power and sample size for linear models of longitudinal data 
        based on the formula due to Liu and Liang (1997) and Diggle et al 
        (1994). Either formula is expressed in terms of marginal model or 
        Generalized Estimating Equations (GEE) parameters.
        This package contains functions which translate pilot mixed 
        effect model parameters (e.g. random intercept and/or slope)
        into marginal model parameters so that the formulas of 
        Diggle et al or Liu and Liang formula can be applied to produce
        sample size calculations for two sample longitudinal designs
        assuming known variance.}
\details{
\tabular{ll}{
Package: \tab longpower\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-06-12\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
%%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Michael C. Donohue <mdonohue@ucsd.edu>
Steven D. Edland
}
\references{
Diggle PJ, Heagerty PJ, Liang K, Zeger SL. Analysis of longitudinal data. Second Edition. Oxford Statistical Science Serires. 2002.

Liu, G., & Liang, K. Y. (1997). Sample size calculations for studies with correlated observations. Biometrics, 53(3), 937-47.
}
\keyword{ package }
\seealso{
\code{lmmpower}
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
%%\examples{
%%~~ simple examples of the most important functions ~~
%%}
