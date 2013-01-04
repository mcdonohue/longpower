diggle.linear.power <-
function(delta, t, sigma2=1, R, 
         sig.level=0.05, power=0.80,
         alternative = c("two.sided", "one.sided"))
{
  alternative = match.arg(alternative)
  if(is.null(dim(R))) R = matrix(R, length(t), length(t)) + diag(1-R,length(t))
  V = sigma2*R
  xi = solve(rbind(1,t)%*%solve(V)%*%cbind(1,t))[2,2]
  n=ceiling(2*(qnorm(ifelse(alternative=="two.sided", sig.level/2, sig.level)) + qnorm(1-power))^2*xi/delta^2)

  NOTE <- "n is number in *each* group"
  METHOD <- "Longitudinal linear model slope power calculation (Diggle et al 2002, page 29)"
  structure(list(n = n, delta = delta, sigma2 = sigma2, R = R, sig.level = sig.level, 
        power = power, alternative = alternative, note = NOTE, 
        method = METHOD), class = "power.longtest")
}

