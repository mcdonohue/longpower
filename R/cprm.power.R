#' Chronic progressive repeated measures (CPRM) model sample size calculations.
#' 
#' This function performs sample size calculations for the chronic progressive 
#' repeated measures (CPRM) model when used to test for differences of change scores 
#' between groups at last visit. Input parameters are random effect variance 
#' and residual error variance as estimated by a REML fit to representative 
#' pilot data or data from a representative prior clinical trial or cohort
#' study.
#' 
#' Default settings perform sample size / power / effect size calculations assuming
#' equal covariance of repeated measures in the 2 groups, equal residual error
#' variance across groups, equal allocation to groups, and assuming no study subject 
#' attrition.  Specifically, variance parameters required for default settings 
#' are `sig2.s`, the variance of random slopes, and `sig2.e`, the residual error
#' variance, both either known or estimated from a mixed model fit by REML
#' to prior data.
#'
#' This function accommodates different variance parameters across groups, 
#' unequal allocation across groups, and study subject attrition (loss to followup), 
#' which may also vary across groups. Details can be found in the description of 
#' \code{\link{edland.linear.power}}
#'
#' @md
#'
#' @name cprm.power
#' @param n sample size, group 1
#' @param lambda allocation ratio (sample size group 1 divided by sample size group 2)
#' @param delta group difference in fixed effect slopes
#' @param t the observation times
#' @param sig2.s variance of random slopes, group 1
#' @param sig2.s_2 variance of random slopes, group 2 (defaults to `sig2.s`)
#' @param sig2.int variance of random intercepts, group 1
#' @param sig2.int_2 variance of random intercepts, group 2 (defaults to `sig2.int`)
#' @param sig.b0b1 covariance of random slopes and intercepts,group 1
#' @param sig.b0b1_2 covariance of random slopes and intercepts, group 2 (defaults to `sig.b0b1`)
#' @param sig2.e residual variance, group 1
#' @param sig2.e_2 residual variance, group 2 (defaults to `sig2.e`)
#' @param p proportion vector for group 1, if i indexes visits, `p[i]` = the proportion whose last visit was at visit `i` (`p` sums to `1`)
#' @param p_2 proportion vector for group 2 (defaults to `p`)
#' @param sig.level type one error
#' @param power power
#' @param alternative one- or two-sided test
#' @param tol	not used (no root finding used in this implementation).
#' @return One of the number of subject required per arm, the `power`, or detectable effect size 
#' given `sig.level` and the other parameter estimates.
#' @author Steven D. Edland, Yu Zhao
#' @seealso [`lmmpower`], [`diggle.linear.power`], [`liu.liang.linear.power`], [`edland.linear.power`], [`hu.mackey.thomas.linear.power`]
#' @references Zhao Y, Edland SD. The chronic progressive repeated measures (CPRM) model for longitudinal data. 
#'  \emph{In process.}
#' @keywords power sample size mixed effects random effects
#' @examples
#' 
#' \dontrun{
#' browseVignettes(package = "longpower")
#' }
#' # An Alzheimer's Disease example using ADAS-cog pilot estimates
#' t <- seq(0,1.5,0.25)
#' cprm.power(delta=1.5, t=t, sig2.s = 24, sig2.e = 10, sig.level=0.05, power = 0.80)
#' 
#' @export cprm.power
#' 
cprm.power <- function(n = NULL, delta = NULL, power = NULL, t = NULL, lambda = 1,
                       sig2.int   = 0,    sig2.s   = NULL, sig.b0b1   = 0,    sig2.e   = NULL, 
                       sig2.int_2 = NULL, sig2.s_2 = NULL, sig.b0b1_2 = NULL, sig2.e_2 = NULL, 
                       sig.level=0.05, p=NULL ,p_2=NULL,
                       alternative = c("two.sided", "one.sided"), tol=NULL) {
  if (sum(sapply(list(n, delta,  power), is.null)) != 1) 
    stop("exactly one of 'n', 'delta', and 'power' must be NULL")
  
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
    stop("'sig.level' must be numeric in [0, 1]")
  
  if (is.null(sig2.s) | is.null(sig2.e)) 
    stop("input values are required for each of sig2.s and sig2.e")
  
  if (is.null(p)) p=rep(c(0,1),c(length(t)-1,1))
  if (is.null(p_2)) p_2=p
  
  if (is.null(sig2.int_2)) sig2.int_2=sig2.int
  if (is.null(sig2.s_2)) sig2.s_2=sig2.s
  if (is.null(sig.b0b1_2)) sig.b0b1_2=sig.b0b1
  if (is.null(sig2.e_2)) sig2.e_2=sig2.e
  
  if (length(t) != length(p) | length(p) != length(p_2))
    stop("t, p, and p_2 must be the same length")
  
  if (sum(p) != 1 | sum(p_2) != 1)
    stop("p and p_2 must sum to 1")
  
  if ( (p[length(p)] !=1 | p_2[length(p)] !=1) & (is.null(sig2.int) | is.null(sig.b0b1)) ) 
    stop("input values are required for each of sig2.int and sig.b0b1 when study subject attrition is specified using p or p_2")
  
  alternative <- match.arg(alternative, c("two.sided", 
                                          "one.sided"))
  
  
  ind_X <- matrix(rep(0,length(t)^2),nrow = length(t))
  ind_X[,1] <- 1
  diag(ind_X)[2:length(t)] <- 1
  
  
  getTerm=function(t,X,sig2.b0,sig.b0b1,sig2.b1,sig2.e,p){
    
    V=matrix(rep(NA,length(t)^2),ncol=length(t))
    for(i in 1:length(t)){
      for(j in 1:length(t)){
        V[i,j]=sig2.int+(t[i]+t[j])*sig.b0b1+t[i]*t[j]*sig2.s 
      }}
    diag(V) = diag(V)+sig2.e
    
    term=0
    k.max=length(t)
    for(k in 2:k.max){
      V.k=V[1:k,1:k]
      X.k=X[1:k,]
      term=term + p[k]*t(X.k)%*%solve(V.k)%*%X.k
    } 
    solve(term)[length(t),length(t)]
  }
  term1=getTerm(t,ind_X,sig2.int  ,sig.b0b1  ,sig2.s  ,sig2.e  ,p)
  term2=getTerm(t,ind_X,sig2.int_2,sig.b0b1_2,sig2.s_2,sig2.e_2,p_2)
  
  
  if (!is.null(n)) 
    n_2 = n/lambda
  
  if (is.null(sig.level) | is.null(sig2.s) | is.null(sig2.e)) 
    stop("solving for sig.level, sig2.s, or sig2.e is no longer supported")
  else if (is.null(power)) {
    power = pnorm(sqrt(n*delta^2/(term1+term2*lambda))+qnorm(ifelse(alternative == "two.sided", sig.level/2, sig.level)))
  } 
  else if (is.null(delta)) {
    delta = sqrt((qnorm(ifelse(alternative == "two.sided", sig.level/2, sig.level))+qnorm(1-power))^2*(term1+term2*lambda)/n)
  } 
  else if (is.null(n)) {
    n = (qnorm(ifelse(alternative == "two.sided", sig.level/2, sig.level))+qnorm(1-power))^2*(term1+term2*lambda)/delta^2
    n_2 = n/lambda
  }
  else # (Shouldn't happen)
    stop("internal error", domain = NA)
  
  METHOD <- "Zhao and Edland, in process"
  structure(list(N = n + n_2, n = c(n, n_2), delta = delta, 
                 t = t, p = p, p_2 = p_2, sig2.int = sig2.int, sig.b0b1 = sig.b0b1, 
                 sig2.s = sig2.s, sig2.e = sig2.e, sig2.int_2 = sig2.int_2, 
                 sig.b0b1_2 = sig.b0b1_2, sig2.s_2 = sig2.s_2, sig2.e_2 = sig2.e_2, 
                 sig.level = sig.level, power = power, alternative = alternative, 
                 note = "N is *total* sample size and n is sample size in *each* group", 
                 method = METHOD), class = "power.longtest")
}
