edland.linear.power <- function(n = NULL, delta = NULL, t = NULL, sig2.s = 0, sig2.e = 1, 
         sig.level=0.05, power=NULL,
         alternative = c("two.sided", "one.sided"))
{
  if (sum(sapply(list(n, delta, sig2.s, sig2.e, power, sig.level), is.null)) != 1) 
      stop("exactly one of 'n', 'delta', 'sig2.s', 'sig2.e', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)

  n.body <- quote({
    ceiling(2*(qnorm(ifelse(alternative=="two.sided", sig.level/2, sig.level)) + qnorm(1-power))^2 * 
      (sig2.s + sig2.e/sum((t - mean(t))^2)) / delta^2)
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
  else if (is.null(sig2.s)) 
      sig2.s <- uniroot(function(sig2.s) eval(n.body) - 
          n, c(1e-10, 1e5))$root
  else if (is.null(sig2.e)) 
      sig2.e <- uniroot(function(sig2.e) eval(n.body) - 
          n, c(1e-10, 1e5))$root
  n <- eval(n.body)

  METHOD <- "Power for longitudinal linear model with random slope (Edland, 2009)"
  structure(list(n = n, delta = delta, sig2.s = sig2.s, sig2.e = sig2.e, sig.level = sig.level, 
        t=t, power = power, alternative = alternative, 
        note = "n is number in *each* group",
        method = METHOD), class = "power.longtest")
}