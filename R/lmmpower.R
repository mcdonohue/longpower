lmmpower <- function(object, ...) UseMethod("lmmpower")
setGeneric("lmmpower")

lmmpower.default <- function(object=NULL,
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
   method = c("diggle", "liuliang"),
   ...)
{
    alternative <- match.arg(alternative)
    method <- match.arg(method)

	n <- length(t)
	if(is.null(R)){
	  D <- matrix(c(sig2.i, cov.s.i, cov.s.i, sig2.s), nrow=2)
	  R = cbind(1,t)%*%D%*%rbind(1,t)  
	  R = R + diag(sig2.e, n, n)
	}
	
	if(method=="liuliang"){
	  u = list(u1 = t, u2 = rep(0,n))
	  v = list(v1 = cbind(1,1,rep(0,n)),
	         v2 = cbind(1,0,t))
	}
	
	if(sum(c(is.null(delta),is.null(pct.change)))!=1) 
		stop("Exactly one of delta and pct.change must be specified.")

	if(is.null(delta)&!is.null(beta))
		delta=pct.change*beta
	else if(is.null(delta)) stop("Parameters beta and pct.change, or delta, must be specified.")
	
	N = ifelse(method=="diggle", 
	      diggle.linear.power(delta, t, R=R, 
	        sig.level=sig.level,
	        power=power,
	        alternative=alternative)$n,
	      liu.liang.linear.power(delta, u=u, v=v, R=R, 
	        sig.level=sig.level,
	        power=power,
	        alternative=alternative)$n)

	if(is.null(delta.lim)&!is.null(beta.lim)) delta.lim <- (delta/beta)*beta.lim
	
	if(!is.null(delta.lim)){
		N.upper = ifelse(method=="diggle",
		            diggle.linear.power(delta.lim[1], t=t, R=R, 
		              sig.level=sig.level,
		              power=power)$n,
		            liu.liang.linear.power(delta.lim[1], u=u, v=v, R=R, 
		              sig.level=sig.level,
		              power=power)$n)
		N.lower = ifelse(method=="diggle",
		            diggle.linear.power(delta.lim[2], t=t, R=R, 
		              sig.level=sig.level,
		              power=power)$n,
		            liu.liang.linear.power(delta.lim[2], u=u, v=v, R=R, 
		              sig.level=sig.level,
		              power=power)$n)
		N.lim <- c(N.lower, N.upper)
		if(N.lim[1]>N.lim[2]) N.lim <- N.lim[2:1]
	}
	else N.lim <- NULL

	results = list(N = N, N.lim = N.lim,
	    beta = beta,
        beta.lim = beta.lim, 
	    pct.change = pct.change, 
        delta = delta, 
        time=t,
        sig2.i=sig2.i, sig2.s = sig2.s, cov.s.i = cov.s.i, sig2.e = sig2.e,
        sig.level = sig.level, power = power, alternative = alternative,
        method = "Power for mixed-effects model",
        R = R
    )
    if(is.character(parameter)){
	    names(results)[3] <- parameter
	    names(results)[4] <- paste(parameter, "lim")
	}
	results <- results[unlist(lapply(results, function(x) !is.null(x)))]
	structure(results, class = "power.longtest")
}

lmmpower.lme <- function(object,
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
{
	
	if(is.numeric(parameter)) parameter <- rownames(summary(object)$tTable)[parameter]
	
	tab <- as.matrix(summary(object$modelStruct)$reStruct[[names(object$groups)]])*object$sigma^2

	if(nrow(tab)>2) stop("Too many random effects. Function is 
	  	equipped to handle at most a random intercept and slope.")
    n = length(t)
	if(is.null(beta))	
	  beta = summary(object)$tTable[parameter,'Value']
	if(is.null(beta.lim))	
	  beta.lim = rep(summary(object)$tTable[parameter,'Value'],2) + 
	    c(-1,1)*qnorm(0.025)*summary(object)$tTable[parameter,'Std.Error']
	if(beta.lim[1] > beta.lim[2]) beta.lim <- beta.lim[2:1]
	# var of random intercept
	if(is.null(sig2.i))
	  sig2.i = as.numeric(tab["(Intercept)","(Intercept)"])
	# var of random slope
	if(is.null(sig2.s))
	  sig2.s = ifelse(nrow(tab)==1, 0, 
	         ifelse(nrow(tab)==2, as.numeric(tab[2,2]), NA))
	# residual var
	if(is.null(sig2.e))
	  sig2.e = object$sigma^2
	# covariance of slope and intercep
	if(is.null(cov.s.i))
	  cov.s.i = ifelse(nrow(tab)==1, 0, 
	          ifelse(nrow(tab)==2, as.numeric(tab[1,2]), NA))

	lmmpower(object=NULL,
		parameter = parameter,
		pct.change = pct.change,
		delta = delta,
		t = t,
		sig.level = sig.level,
		power = power, 
		alternative = alternative,
		beta=beta,
		beta.lim=beta.lim,
		delta.lim=delta.lim,
		sig2.i=sig2.i,
		sig2.s=sig2.s,
		sig2.e=sig2.e,
		cov.s.i=cov.s.i, ...)
}

lmmpower.gee <- function(object,
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
{
	
	if(is.numeric(parameter)) parameter <- rownames(summary(object)$coefficients)[parameter]
	
	if(is.null(beta))	
	  beta = summary(object)$coefficients[parameter,'Estimate']
	if(is.null(beta.lim))	
	  beta.lim = rep(summary(object)$coefficients[parameter,'Estimate'],2) + 
	    c(-1,1)*qnorm(0.025)*summary(object)$coefficients[parameter,'Robust S.E.']
	if(beta.lim[1] > beta.lim[2]) beta.lim <- beta.lim[2:1]
	
	R <- summary(object)$working.correlation * summary(object)$scale

	lmmpower(object=NULL,
		parameter = parameter,
		pct.change = pct.change,
		sig.level = sig.level,
		power = power, 
		alternative = alternative,
		beta=beta,
		beta.lim=beta.lim,
		delta.lim=delta.lim,
		R=R,
		t=t, ...)
}

setMethod("lmmpower", signature(object = "mer"),
  function(object,
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
{
	if(is.numeric(parameter)) parameter <- rownames(summary(object)@coefs)[parameter]

	tab <- summary(object)@REmat
	if(nrow(tab)>3) stop("Too many random effects. Function is 
	  	equipped to handle at most a random intercept and slope.")
    n = length(t)
	if(is.null(beta))	
	  beta = coef(summary(object))[parameter,'Estimate']
	if(is.null(beta.lim))	
	  beta.lim = rep(coef(summary(object))[parameter,"Estimate"],2) + 
	    c(-1,1)*qnorm(0.025)*coef(summary(object))[parameter,"Std. Error"]
	if(beta.lim[1] > beta.lim[2]) beta.lim <- beta.lim[2:1]
	# var of random intercept
	if(is.null(sig2.i))
	  sig2.i = as.numeric(tab[1,"Variance"])
	# var of random slope
	if(is.null(sig2.s))
	  sig2.s = ifelse(nrow(tab)==2, 0, 
	         ifelse(nrow(tab)==3, as.numeric(tab[2,"Variance"]), NA))
	# residual var
	if(is.null(sig2.e))
	  sig2.e = as.numeric(tab[tab[,"Groups"]=="Residual","Variance"])
	# covariance of slope and intercep
	if(is.null(cov.s.i))
	  cov.s.i = ifelse(nrow(tab)==2, 0, 
	          ifelse(nrow(tab)==3, as.numeric(tab[2,"Corr"])*sqrt(sig2.i)*sqrt(sig2.s), NA))
	lmmpower(object=NULL,
		parameter = parameter,
		pct.change = pct.change,
		delta = delta,
		t = t,
		sig.level = sig.level,
		power = power, 
		alternative = alternative,
		beta=beta,
		beta.lim=beta.lim,
		delta.lim=delta.lim,
		sig2.i=sig2.i,
		sig2.s=sig2.s,
		sig2.e=sig2.e,
		cov.s.i=cov.s.i, ...)
})
