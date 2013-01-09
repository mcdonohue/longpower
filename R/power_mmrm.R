power.mmrm <- function(Ra, ra, sigmaa, Rb = NULL, rb = NULL, sigmab = NULL, lambda = 1,
  delta = NULL, sig.level = 0.05, power = NULL, 
  alternative = c("two.sided", "one.sided"))
{
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)

  # formula (3) on page 4 in
  # Lu, K., Luo, X., & Chen, P.-Y. (July 14, 2008). Sample size estimation for repeated measures analysis in randomized clinical trials with missing data. International Journal of Biostatistics, 4, (1)
    
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
    phib <- phia
  }
  if(is.null(sigmab)){
    sigma <- sigmaa
  }else{
    sigma <- mean(c(sigmaa, sigmab))
  }
  
  Na <- as.numeric(ceiling(
    (phia + lambda * phib)*(qnorm(ifelse(alternative=="two.sided", sig.level/2, sig.level)) + qnorm(1-power))^2*
    sigma^2/delta^2
  ))
  Nb <- as.numeric(ceiling(Na/lambda))
  
  METHOD <- "Power for Mixed Model of Repeated Measures (Lu, Luo, & Chen, 2008)"
  structure(list(n1 = Na, n2 = Nb, delta = delta, sig.level = sig.level, 
        power = power, alternative = alternative, 
        method = METHOD), class = "power.htest")
}

power.mmrm.ar1 <- function(rho, ra, sigmaa, rb = NULL, sigmab = NULL, 
  lambda = 1, times = 1:length(ra),
  delta = NULL, sig.level = 0.05, power = NULL, 
  alternative = c("two.sided", "one.sided"))
{
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)

  # formula (3) on page 4 in
  # Lu, K., Luo, X., & Chen, P.-Y. (July 14, 2008). Sample size estimation for repeated measures analysis in randomized clinical trials with missing data. International Journal of Biostatistics, 4, (1)
  if(length(times) != length(ra)) stop("ra and times should be the same length")
  
  J <- length(ra)
  phia <- 1/ra[J] - sum(
    rho^(2*(times[J] - times[-J])) *
    (1/ra[-1] - 1/ra[-J]))

  if(!is.null(rb)){
    phib <- 1/rb[J] - sum(
      rho^(2*(times[J] - times[-J])) *
      (1/rb[-1] - 1/rb[-J]))
  }else{
    phib <- phia
  }
  if(is.null(sigmab)){
    sigma <- sigmaa
  }else{
    sigma <- mean(c(sigmaa, sigmab))
  }
  
  Na <- as.numeric(ceiling(
    (phia + lambda * phib)*(qnorm(ifelse(alternative=="two.sided", sig.level/2, sig.level)) + qnorm(1-power))^2*
    sigma^2/delta^2
  ))
  Nb <- as.numeric(ceiling(Na/lambda))
  
  METHOD <- "Power for Mixed Model of Repeated Measures (Lu, Luo, & Chen, 2008)"
  structure(list(n1 = Na, n2 = Nb, phi1 = phia, phi2 = phib, 
        delta = delta, times = times, sig.level = sig.level, 
        power = power, alternative = alternative, 
        method = METHOD), class = "power.htest")
}
