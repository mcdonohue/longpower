diggle.linear.power <-
function(n=NULL, delta=NULL, t=NULL, sigma2=1, R=NULL, 
         sig.level=0.05, power=NULL,
         alternative=c("two.sided", "one.sided"))
{
  if (sum(sapply(list(n, delta, sigma2, power, sig.level), is.null)) != 1) 
      stop("exactly one of 'delta', 'sigma2', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)

  if(is.null(dim(R))) R = matrix(R, length(t), length(t)) + diag(1-R,length(t))
  
  n.body <- quote({
    V = sigma2*R
    xi = solve(rbind(1,t)%*%solve(V)%*%cbind(1,t))[2,2]
    2*(qnorm(ifelse(alternative=="two.sided", sig.level/2, sig.level)) +
       qnorm(1-power))^2*xi/delta^2
  })
  
  if (is.null(sig.level)) 
      sig.level <- uniroot(function(sig.level) eval(n.body) - 
          n, c(1e-10, 1 - 1e-10))$root
  else if (is.null(power)) 
      power <- uniroot(function(power) eval(n.body) - 
          n, c(1e-3, 1 - 1e-10))$root
  else if (is.null(delta)) 
      delta <- uniroot(function(delta) eval(n.body) - 
          n, c(1e-10, 1e5))$root
  else if (is.null(sigma2)) 
      sigma2 <- uniroot(function(sigma2) eval(n.body) - 
          n, c(1e-10, 1e5))$root
  n <- eval(n.body)

  METHOD <- "Longitudinal linear model slope power calculation (Diggle et al 2002, page 29)"
  structure(list(n = n, delta = delta, sigma2 = sigma2, R = R, sig.level = sig.level, 
        power = power, alternative = alternative, 
        note = "n is number in *each* group",
        method = METHOD), class = "power.longtest")
}
