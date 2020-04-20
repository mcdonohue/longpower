
#' @export
lmmpower <- function(object, ...) UseMethod("lmmpower")
setGeneric("lmmpower")

#' Sample size calculations for linear mixed models of rate of change based on
#' lmer, lme, or gee "placebo" pilot estimates.
#' 
#' These functions compute sample size for linear mixed models based on the
#' formula due to Diggle (2002) or Liu and Liang (1997).  These formulae are
#' expressed in terms of marginal model or Generalized Estimating Equations
#' (GEE) parameters.  These functions translate pilot mixed effect model
#' parameters (e.g. random intercept and/or slope, fixed effects, etc.)  into
#' marginal model parameters so that either formula can be applied to
#' equivalent affect. Pilot estimates are assumed to be from an appropriate
#' "placebo" group and the parameter of interest is assumed to be the rate of
#' change over time of the outcome.
#' 
#' Any parameters not explicitly stated are extracted from the fitted
#' \code{object}.
#' 
#' @name lmmpower
#' @aliases lmmpower-methods lmmpower,ANY-method lmmpower,merMod-method
#' lmmpower.default lmmpower.lme lmmpower.gee lmmpower.numeric lmmpower.double
#' @docType methods
#' @param object an object returned by lme4
#' @param n sample size per group
#' of a mixed-effects model object to placebo data assumed to have either a
#' random intercept, or a random intercept and random effect for time (slope);
#' and fixed effect representing the rate of change in a placebo group.
#' @param parameter the name or position
#' of the rate of change parameter of interest, e.g. (\code{"time"},
#' \code{"t"}, or \code{2} if it is the second specified fixed effect).
#' @param pct.change the percent change
#' in the pilot estimate of the parameter of interest (\code{beta}, the
#' placebo/null effect)
#' @param delta the change in the pilot estimate
#' of the parameter of interest, computed from \code{pct.change} if left
#' missing.
#' @param t vector of time points
#' @param sig.level Type I error
#' @param power power
#' @param alternative \code{"two.sided"} or \code{"one.sided"}
#' @param beta pilot estimate of the placebo
#' effect (slope or rate of change in the outcome)
#' @param beta.CI 95\% confidence limits of
#' the pilot estimate of beta
#' @param delta.CI 95\% confidence limits of
#' the effect size
#' @param sig2.i pilot estimate of variance
#' of random intercept
#' @param sig2.s pilot estimate of variance
#' of random slope
#' @param sig2.e pilot estimate of residual
#' variance
#' @param cov.s.i pilot estimate of
#' covariance of random slope and intercept
#' @param R pilot estimate of a marginal
#' model working correlation matrix
#' @param method the formula to use. Defaults
#' to \code{"diggle"} for Diggle et al (2002). Alternatively \code{"liuliang"}
#' can be selected for Liu & Liang (1997), or  \code{"edland"} for Ard & Edland (2011).
#' @param tol numerical tolerance used in root finding.
#' @param ... other arguments
#' @return An object of class \code{power.htest} giving the calculated sample
#' size, N, per group and other parameters.
#' @author Michael C. Donohue
#' @seealso \code{\link{liu.liang.linear.power}},
#' \code{\link{diggle.linear.power}}, \code{\link{edland.linear.power}}
#' @references Diggle P.J., Heagerty P.J., Liang K., Zeger S.L. (2002)
#' \emph{Analysis of longitudinal data}. Second Edition. Oxford Statistical
#' Science Series.
#' 
#' Liu, G., and Liang, K. Y. (1997) Sample size calculations for studies with
#' correlated observations. \emph{Biometrics}, 53(3), 937-47.
#' 
#' Ard, C. and Edland, S.D. (2011) Power calculations for clinical trials in Alzheimer's disease. 
#' \emph{Journal of Alzheimer's Disease.} 21:369-377. 
#' 
#' @keywords power sample size mixed effects random effects marginal model
#' methods
#' @examples
#' 
#' \dontrun{
#' browseVignettes(package = "longpower")
#' }
#' 
#' lmmpower(delta=1.5, t = seq(0,1.5,0.25),
#' 	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)
#' lmmpower(n=208, t = seq(0,1.5,0.25),
#' 	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)
#' lmmpower(beta = 5, pct.change = 0.30, t = seq(0,1.5,0.25),
#' 	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)
#' 
#' \dontrun{
#' library(lme4)
#' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
#' lmmpower(fm1, pct.change = 0.30, t = seq(0,9,1), power = 0.80)
#' 
#' library(nlme)
#' fm2 <- lme(Reaction ~ Days, random=~Days|Subject, sleepstudy)
#' lmmpower(fm2, pct.change = 0.30, t = seq(0,9,1), power = 0.80)
#' 
#' # random intercept only
#' fm3 <- lme(Reaction ~ Days, random=~1|Subject, sleepstudy)
#' lmmpower(fm3, pct.change = 0.30, t = seq(0,9,1), power = 0.80)
#' 
#' library(gee)
#' fm4 <- gee(Reaction ~ Days, id = Subject,
#'             data = sleepstudy,
#'             corstr = "exchangeable")
#' lmmpower(fm4, pct.change = 0.30, t = seq(0,9,1), power = 0.80)
#' }
#' 
#' @method lmmpower default
#' @export
lmmpower.default <- function(object=NULL,
   n=NULL,
   parameter = 2,
   pct.change = NULL,
   delta = NULL,
   t = NULL,
   sig.level = 0.05,
   power = NULL, 
   alternative = c("two.sided", "one.sided"),
   beta=NULL,
   beta.CI=NULL,
   delta.CI=NULL,
   sig2.i=NULL,
   sig2.s=NULL,
   sig2.e=NULL,
   cov.s.i=NULL,
   R=NULL,
   method = c("diggle", "liuliang", "edland"),
   tol = .Machine$double.eps^2,
   ...)
{
  if(sum(!sapply(list(delta, pct.change), is.null))==2) 	
		stop("Only one of delta and pct.change must be specified.")
	if(is.null(delta)&!is.null(beta)&!is.null(pct.change))
	  delta<-pct.change*beta
  if (sum(sapply(list(n, delta, power, sig.level), is.null)) != 1) 
      stop("exactly one of 'n', 'delta', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  
	m <- length(t)
	if(is.null(R) & method %in% c("diggle", "liuliang")){
	  D <- matrix(c(sig2.i, cov.s.i, cov.s.i, sig2.s), nrow=2)
	  R <- cbind(1,t)%*%D%*%rbind(1,t)  
	  R <- R + diag(sig2.e, m, m)
	}

	if(method == "liuliang"){
	  u <- list(u1 = t, u2 = rep(0,m))
	  v <- list(v1 = cbind(1,1,rep(0,m)),
	         v2 = cbind(1,0,t))
	  if(!is.null(n)) N <- n*2 else N <- NULL
	}

	results <- switch(method,
	  edland = edland.linear.power(n=n, delta=delta, t=t, 
      sig2.s=sig2.s, sig2.e=sig2.e, 
      sig.level=sig.level,
      power=power,
      alternative=alternative,tol=tol,...),
	  diggle = diggle.linear.power(n=n, delta=delta, t=t, R=R, 
	    sig.level=sig.level,
	    power=power,
	    alternative=alternative,tol=tol,...),
	  liuliang = liu.liang.linear.power(N=N, delta=delta, u=u, v=v, R=R,
	    sig.level=sig.level,
	    power=power,
	    alternative=alternative,tol=tol,...))

	if(is.null(delta.CI)&!is.null(beta.CI)) results$delta.CI <- (results$delta/beta)*beta.CI
	if(!is.null(beta)) results$beta <- beta
	if(!is.null(beta.CI)) results$beta.CI <- beta.CI
		
	if(!is.null(results$delta.CI)){
		n.upper <- switch(method,
		  edland = edland.linear.power(n=NULL, results$delta.CI[1], t=t, 
		    sig2.s=sig2.s, sig2.e=sig2.e, 
		    sig.level=sig.level,
		    power=power,
		    alternative=alternative,tol=tol,...)$n,
      diggle = diggle.linear.power(n=NULL, results$delta.CI[1], t=t, R=R, 
        sig.level=sig.level,
        power=power,
        alternative=alternative,tol=tol,...)$n,
		  liuliang = liu.liang.linear.power(N=NULL, results$delta.CI[1], u=u, v=v, R=R, 
		    sig.level=sig.level,
		    power=power,tol=tol,...)$N/2)
		n.lower <- switch(method,
		  edland = edland.linear.power(n=NULL, results$delta.CI[2], t=t, 
		    sig2.s=sig2.s, sig2.e=sig2.e, 
        sig.level=sig.level,
        power=power,
        alternative=alternative,tol=tol,...)$n,
      diggle = diggle.linear.power(n=NULL, results$delta.CI[2], t=t, R=R, 
		    sig.level=sig.level,
		    power=power,tol=tol,...)$n,
		  liuliang = liu.liang.linear.power(N=NULL, results$delta.CI[2], u=u, v=v, R=R, 
		    sig.level=sig.level,
		    power=power,tol=tol,...)$N/2)
		n.CI <- c(n.lower, n.upper)
		if(n.CI[1]>n.CI[2]) n.CI <- n.CI[2:1]
		results$n.CI <- n.CI 
	}

    if(is.character(parameter)){
      names(results)[names(results) == "beta"] <- parameter
      names(results)[names(results) == "beta.CI"] <- paste(parameter, "CI")
  }
	results <- results[unlist(lapply(results, function(x) !is.null(x)))]
	structure(results, class = "power.longtest")
}

#' @export
#' @method lmmpower double
lmmpower.double <- lmmpower.default

#' @export
#' @method lmmpower numeric
lmmpower.numeric <- lmmpower.default

#' @importFrom nlme getVarCov
#' @method lmmpower lme
#' @export
lmmpower.lme <- function(object,
   n = NULL,
   parameter = 2,
   pct.change = NULL,
   delta = NULL,
   t = NULL,
   sig.level = 0.05,
   power = NULL, 
   alternative = c("two.sided", "one.sided"),
   beta=NULL,
   beta.CI=NULL,
   delta.CI=NULL,
   sig2.i=NULL,
   sig2.s=NULL,
   sig2.e=NULL,
   cov.s.i=NULL,
   method = c("diggle", "liuliang", "edland"),
   tol = .Machine$double.eps^2,
   ...)
{
	alternative <- match.arg(alternative)
  method <- match.arg(method)
  
	if(is.numeric(parameter)) parameter <- rownames(summary(object)$tTable)[parameter]
	
	tab <- getVarCov(object)

	if(nrow(tab)>2) stop("Too many random effects. Function is 
	  	equipped to handle at most a random intercept and slope.")

	if(is.null(beta))	
	  beta = summary(object)$tTable[parameter,'Value']
	if(is.null(beta.CI))	
	  beta.CI = rep(summary(object)$tTable[parameter,'Value'],2) + 
	    c(-1,1)*qnorm(0.025)*summary(object)$tTable[parameter,'Std.Error']
	if(beta.CI[1] > beta.CI[2]) beta.CI <- beta.CI[2:1]
	# var of random intercept
	if(is.null(sig2.i))
	  sig2.i = tab["(Intercept)", "(Intercept)"]
	# var of random slope
	if(is.null(sig2.s))
	  sig2.s = ifelse(nrow(tab)==1, 0, 
	         ifelse(nrow(tab)==2, tab[2, 2], NA))
	
	# residual var
	if(is.null(sig2.e))
	  sig2.e = object$sigma^2
	# covariance of slope and intercep
	if(is.null(cov.s.i))
	  cov.s.i = ifelse(nrow(tab)==1, 0, 
	          ifelse(nrow(tab)==2, tab[1, 2], NA))

	lmmpower(object=NULL,
		n = n,
		parameter = parameter,
		pct.change = pct.change,
		delta = delta,
		t = t,
		sig.level = sig.level,
		power = power, 
		alternative = alternative,
		beta = beta,
		beta.CI = beta.CI,
		delta.CI = delta.CI,
		sig2.i = sig2.i,
		sig2.s = sig2.s,
		sig2.e = sig2.e,
		cov.s.i = cov.s.i, 
		method = method,
    tol=tol, ...)
}

#' @method lmmpower gee
#' @export
lmmpower.gee <- function(object,
   n = NULL,
   parameter = 2,
   pct.change = NULL,
   delta = NULL,
   t = NULL,
   sig.level = 0.05,
   power = NULL, 
   alternative = c("two.sided", "one.sided"),
   beta=NULL,
   beta.CI=NULL,
   delta.CI=NULL,
   method = c("diggle", "liuliang", "edland"),
   tol = .Machine$double.eps^2,
   ...)
{
	alternative <- match.arg(alternative)
  method <- match.arg(method)
  
	if(is.numeric(parameter)) parameter <- rownames(summary(object)$coefficients)[parameter]
	
	if(is.null(beta))	
	  beta = summary(object)$coefficients[parameter,'Estimate']
	if(is.null(beta.CI))	
	  beta.CI = rep(summary(object)$coefficients[parameter,'Estimate'],2) + 
	    c(-1,1)*qnorm(0.025)*summary(object)$coefficients[parameter,'Robust S.E.']
	if(beta.CI[1] > beta.CI[2]) beta.CI <- beta.CI[2:1]
	
	R <- summary(object)$working.correlation * summary(object)$scale

	lmmpower(object=NULL,
		n = n,
		parameter = parameter,
		pct.change = pct.change,
		sig.level = sig.level,
		power = power, 
		alternative = alternative,
		sig2.e=1,
		sig2.s=0,
		beta=beta,
		beta.CI=beta.CI,
		delta.CI=delta.CI,
		R=R,
		t=t, 
		method=method,
    tol=tol, ...)
}

#' @importFrom lme4 VarCorr fixef getME
#' @export
setMethod("lmmpower", signature(object = "merMod"),
  function(object, 
   n = NULL, 
   parameter = 2,
   pct.change = NULL,
   delta = NULL,
   t = NULL,
   sig.level = 0.05,
   power = NULL, 
   alternative = c("two.sided", "one.sided"),
   beta=NULL,
   beta.CI=NULL,
   delta.CI=NULL,
   sig2.i=NULL,
   sig2.s=NULL,
   sig2.e=NULL,
   cov.s.i=NULL,
   method = c("diggle", "liuliang", "edland"),
   tol = .Machine$double.eps^2,
   ...)
{
  if (!(sum(sapply(list(n, delta, power, sig.level), is.null)) == 1 |
        sum(sapply(list(n, pct.change, power, sig.level), is.null)) == 1)) 
      stop("exactly one of 'n', 'delta' (or 'pct.change'), 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
	alternative <- match.arg(alternative)
  method <- match.arg(method)
  
	if(is.numeric(parameter)) parameter <- rownames(coef(summary(object)))[parameter]
	
	tab <- VarCorr(object)
	if(length(tab)>1) stop("Too many grouping levels. Function is 
	  	equipped to handle at most one grouping level.")
	tab <- tab[[1]]

	if(nrow(tab)>3) stop("Too many random effects. Function is 
	  	equipped to handle at most a random intercept and slope.")
  m <- length(t)
	if(is.null(beta))	
	  beta = fixef(object)[parameter][[1]]
	if(is.null(beta.CI))	
	  beta.CI = rep(coef(summary(object))[parameter,"Estimate"],2) + 
	    c(-1,1)*qnorm(0.025)*coef(summary(object))[parameter,"Std. Error"]
	if(beta.CI[1] > beta.CI[2]) beta.CI <- beta.CI[2:1]

	# var of random intercept
	if(is.null(sig2.i))
	  sig2.i = tab[1, 1]
	# var of random slope
	if(is.null(sig2.s))
	  sig2.s = ifelse(nrow(tab)==1, 0, 
	         ifelse(nrow(tab)==2, tab[2 ,2], NA))
	# residual var
	if(is.null(sig2.e))
	  sig2.e = getME(object, "sigma")^2
	# covariance of slope and intercept
	if(is.null(cov.s.i))
    cov.s.i = ifelse(nrow(tab)==1, 0, 
            ifelse(nrow(tab)==2, tab[2, 1], NA))

	lmmpower(n=n, object=NULL,
		parameter = parameter,
		pct.change = pct.change,
		delta = delta,
		t = t,
		sig.level = sig.level,
		power = power, 
		alternative = alternative,
		beta=beta,
		beta.CI=beta.CI,
		delta.CI=delta.CI,
		sig2.i=sig2.i,
		sig2.s=sig2.s,
		sig2.e=sig2.e,
		cov.s.i=cov.s.i, 
		method = method, 
    tol=tol, ...)
})
