liu.liang.linear.power <- function(n=NULL, delta=NULL, u=NULL, v=NULL, sigma2=1, R=NULL, R.list=NULL,
  sig.level=0.05, power=NULL, 
  Pi = rep(1/length(u),length(u)),
  alternative = c("two.sided", "one.sided"))
{
  if (sum(sapply(list(n, delta, sigma2, power, sig.level), is.null)) != 1) 
      stop("exactly one of 'n', 'sigma2', 'delta', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)

  if(sum(c(!is.null(R), !is.null(R.list))) != 1) 
    stop("Exactly one of R or R.list must be specified.")
  if(sum(Pi) != 1) 
    stop("Pi must sum to 1.")

  if(!is.null(R)){ # make list of Rs
    R.list <- lapply(1:length(u), function(i) R)
  }

  # invert R.list
  Rinv <- lapply(1:length(R.list), 
    function(i){
      R <- R.list[[i]]
      # if R is not a matrix, we assume exchangeable correlation structure
      if(is.null(dim(R)) & length(R) == 1 & length(u[[i]]) > 1){
        R <- matrix(R, length(u[[i]]), length(u[[i]])) + diag(1-R,length(u[[i]]))
      }else if(is.null(dim(R)) & length(R) == 1 & length(u[[i]]) == 1){
        R <- matrix(R, length(u[[i]]), length(u[[i]]))
      }
      return(solve(R))
    }
  )

  n.body <- quote({
    Ipl <- 0
    for(i in 1:length(u)) 
    Ipl <- Ipl + Pi[i]*t(u[[i]])%*%Rinv[[i]]%*%v[[i]]
    Ipl <- Ipl/sigma2   

    Ill <- 0
    for(i in 1:length(u)) 
    Ill <- Ill + Pi[i]*t(v[[i]])%*%Rinv[[i]]%*%v[[i]]
    Illinv <- solve(Ill/sigma2)

    Sigma1 <- 0
    for(i in 1:length(u)) 
    Sigma1 <- Sigma1 + Pi[i]*(t(u[[i]])-Ipl%*%Illinv%*%t(v[[i]]))%*%Rinv[[i]]%*%
    (u[[i]]-v[[i]]%*%Illinv%*%t(Ipl))                           
    Sigma1 <- Sigma1/sigma2

    ceiling((qnorm(1-ifelse(alternative=="two.sided", sig.level/2, sig.level)) + 
    qnorm(power))^2/
    (delta%*%Sigma1%*%delta)[1,1])
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

  METHOD <- "Longitudinal linear model power calculation (Liu & Liang, 1997)"
  structure(list(n = n, delta = delta, sigma2 = sigma2, 
    sig.level = sig.level, power = power, alternative = alternative,
    R = R,
    note = "n is *1/2* the total sample size required.",
    method = METHOD), class = "power.longtest")
}