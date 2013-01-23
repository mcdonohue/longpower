power.mmrm <- function(N = NULL, Ra = NULL, ra = NULL, sigmaa = NULL, 
  Rb = NULL, rb = NULL, sigmab = NULL, lambda = 1,
  delta = NULL, sig.level = 0.05, power = NULL, 
  alternative = c("two.sided", "one.sided"))
{
  if (sum(sapply(list(N, delta, power, sig.level), is.null)) != 1) 
      stop("exactly one of 'N', 'delta', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)

  # formula (3) on page 4 in
  # Lu, K., Luo, X., & Chen, P.-Y. (July 14, 2008). Sample size estimation for repeated measures analysis in randomized clinical trials with missing data. International Journal of Biostatistics, 4, (1)
    
  n.body <- quote({
    Ia <- 0
    ra <- c(ra, 0)
    for(j in 1:nrow(Ra)){
      Raj <- matrix(0, nrow(Ra), nrow(Ra))
      Raj[1:j, 1:j] <- solve(Ra[1:j, 1:j])
      Ia <- Ia + (ra[j] - ra[j+1]) * Raj
    }
    phia <- solve(Ia)[j,j]
    if(!is.null(Rb) & !is.null(rb)){
      Ib <- 0
      rb <- c(rb, 0)
      for(j in 1:nrow(Rb)){
        Rbj <- matrix(0, nrow(Rb), nrow(Rb))
        Rbj[1:j, 1:j] <- solve(Rb[1:j, 1:j])
        Ib <- Ib + (rb[j] - rb[j+1]) * Rbj
      }
      phib <- solve(Ib)[j,j]  
    }else{
      rb <- ra
      phib <- phia
    }
    if(is.null(sigmab)){
      sigma <- sigmaa
    }else{
      sigma <- mean(c(sigmaa, sigmab))
    }
  
    Na <- as.numeric(
      (phia + lambda * phib)*(qnorm(ifelse(alternative=="two.sided", sig.level/2, sig.level)) + qnorm(1-power))^2*
      sigma^2/delta^2
    )
    Nb <- as.numeric(Na/lambda)
    list(N = Na + Nb, Na = Na, Nb = Nb)
  })
    
  if (is.null(sig.level)) 
      sig.level <- uniroot(function(sig.level) eval(n.body)$N - 
          N, c(1e-10, 1 - 1e-10))$root
  else if (is.null(power)) 
      power <- uniroot(function(power) eval(n.body)$N - 
          N, c(1e-3, 1 - 1e-10))$root
  else if (is.null(delta)) 
      delta <- uniroot(function(delta) eval(n.body)$N - 
          N, c(1e-10, 1e5))$root
    
  Ns <- eval(n.body)
  Na <- Ns$Na
  Nb <- Ns$Nb
  METHOD <- "Power for Mixed Model of Repeated Measures (Lu, Luo, & Chen, 2008)"
  structure(list(n1 = Na, n2 = Nb, 
        retention1 = ra[-length(ra)], retention2 = rb[-length(rb)], 
        delta = delta, sig.level = sig.level, 
        power = power, alternative = alternative, 
        method = METHOD), class = "power.htest")
}

power.mmrm.ar1 <- function(N = NULL, rho = NULL, 
  ra = NULL, sigmaa = NULL, rb = NULL, sigmab = NULL, 
  lambda = 1, times = 1:length(ra),
  delta = NULL, sig.level = 0.05, power = NULL, 
  alternative = c("two.sided", "one.sided"))
{
  if (sum(sapply(list(N, rho, delta, power, sig.level), is.null)) != 1) 
      stop("exactly one of 'N', 'rho', 'delta', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)

  # formula (3) on page 4 in
  # Lu, K., Luo, X., & Chen, P.-Y. (July 14, 2008). Sample size estimation for repeated measures analysis in randomized clinical trials with missing data. International Journal of Biostatistics, 4, (1)
  if(length(times) != length(ra)) stop("ra and times should be the same length")
  
  phia <- phib <- NULL
  
  n.body <- quote({
      J <- length(ra)
      phia <- 1/ra[J] - sum(
        rho^(2*(times[J] - times[-J])) *
        (1/ra[-1] - 1/ra[-J]))

      if(!is.null(rb)){
        phib <- 1/rb[J] - sum(
          rho^(2*(times[J] - times[-J])) *
          (1/rb[-1] - 1/rb[-J]))
      }else{
        rb <- ra
        phib <- phia
      }
      if(is.null(sigmab)){
        sigma <- sigmaa
      }else{
        sigma <- mean(c(sigmaa, sigmab))
      }
  
      Na <- as.numeric(
        (phia + lambda * phib)*(qnorm(ifelse(alternative=="two.sided", sig.level/2, sig.level)) + qnorm(1-power))^2*
        sigma^2/delta^2
      )
      Nb <- as.numeric(Na/lambda)
      list(N = Na + Nb, Na = Na, Nb = Nb)
  })
  
  if (is.null(sig.level)) 
      sig.level <- uniroot(function(sig.level) eval(n.body)$N - 
          N, c(1e-10, 1 - 1e-10))$root
  else if (is.null(power)) 
      power <- uniroot(function(power) eval(n.body)$N - 
          N, c(1e-3, 1 - 1e-10))$root
  else if (is.null(delta)) 
      delta <- uniroot(function(delta) eval(n.body)$N - 
          N, c(1e-10, 1e5))$root
  else if (is.null(rho)) 
      rho <- uniroot(function(rho) eval(n.body)$N - 
          N, c(1e-10, 1 - 1e-10))$root
    
  Ns <- eval(n.body)
  Na <- Ns$Na
  Nb <- Ns$Nb
  
  METHOD <- "Power for Mixed Model of Repeated Measures (Lu, Luo, & Chen, 2008)"
  structure(list(n1 = Na, n2 = Nb, rho = rho, 
        retention1 = ra[-length(ra)], retention2 = rb[-length(rb)],
        phi1 = phia, phi2 = phib, 
        delta = delta, times = times, sig.level = sig.level, 
        power = power, alternative = alternative, 
        method = METHOD), class = "power.htest")
}