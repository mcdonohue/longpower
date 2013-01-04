% Adapted form Peter Dalgaard's print.power.htest

\name{print.power.longtest}
\alias{print.power.longtest}
\title{Print method for longitudinal data power calculation object}
\description{
  Print object of class \code{"power.longtest"} in nice layout.
}
\usage{
\method{print}{power.longtest}(x, \dots)
}
\arguments{
  \item{x}{Object of class \code{"power.longtest"}.}
  \item{\dots}{further arguments to be passed to or from methods.}  
}
\details{
  A \code{power.longtest} object is just a named list of numbers and
  character strings, supplemented with \code{method} and \code{note}
  elements.  The \code{method} is displayed as a title, the \code{note}
  as a footnote, and the remaining elements are given in an aligned
  \sQuote{name = value} format.
}
\value{
  none
}

\seealso{
  \code{\link{liu.liang.linear.power}},
  \code{\link{diggle.linear.power}},
  \code{\link{lmmpower}},
}
\keyword{longtest}
