#' Power for two-stage mixed effects model
#'
#' The two-stage mixed effects model details can be found
#' in Section 8.4 of Fitzmaurice, Laird, and Ware (2011).
#'
#' This program was co-developed by Nan Hu (Genentech Inc., Biostatistics)
#' and a summer intern Zhe Qu (Ph.D. candidate at Tulane University). If you
#' have any technical questions or any ideas for improving functionality,
#' please contact Nan Hu (hu.nan@@gene.com).
#'
#' @param n sample size, group 1
#' @param lambda allocation ratio (= (sample size group 1)/(sample size group 2)
#' @param delta Effect size (absolute difference in rate of decline between tx and placebo)
#' @param duration Duration of follow-up
#' @param n_assessments Number of assessments (including baseline)
#' @param sig2.b Between-subject variance
#' @param sig2.w Within-subject variance
#' @param sig.level type one error
#' @param power power
#' @param alternative one- or two-sided test
#' @param tol	numerical tolerance used in root finding
#' @return One of the number of subject required per arm, the `power`, or
#' detectable effect size given `sig.level` and the other parameter estimates.
#'
#' @author Nan Hu, Michael C. Donohue
#' @seealso \code{\link{lmmpower}}, \code{\link{edland.linear.power}}, 
#' \code{\link{rcrm.power}}
#' @references Hu, N., Mackey, H., & Thomas, R. (2021). Power and sample size
#' for random coefficient regression models in randomized experiments with
#' monotone missing data. \emph{Biometrical Journal}, 63(4), 806-824.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' browseVignettes(package = "longpower")
#' }
#' # An Alzheimer's Disease example using ADAS-cog pilot estimates
#' t <- seq(0,1.5,0.25)
#'
#' two.stage.me.power(delta=1.5, duration = 1.5, n_assessments = 7, sig2.b = 24, sig2.w = 10, sig.level=0.05, power = 0.80)
#' two.stage.me.power(n=207, duration = 1.5, n_assessments = 7, sig2.b = 24, sig2.w = 10, sig.level=0.05, power = 0.80)
#' two.stage.me.power(n=207, delta=1.5, duration = 1.5, n_assessments = 7, sig2.b = 24, sig2.w = 10, sig.level=0.05)
#'
two.stage.me.power <-
  function(n = NULL,
    delta = NULL,
    power = NULL,
    duration = NULL,
    n_assessments = NULL,
    lambda = 1,
    sig2.b = 1,
    sig2.w = 1,
    sig.level = 0.05,
    alternative = c("two.sided", "one.sided"),
    tol = .Machine$double.eps^2)
  {
    # adapted from https://github.com/nan-hu-personal/RCRM_power_size/blob/master/app.R
  
    if (sum(sapply(list(n, delta,  power), is.null)) != 1)
      stop("exactly one of 'n', 'delta', and 'power' must be NULL")
    
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 >
        sig.level | sig.level > 1))
      stop("'sig.level' must be numeric in [0, 1]")
    
    if (!is.null(power) && !is.numeric(power) || any(0 >
        power | power > 1))
      stop("'power' must be numeric in [0, 1]")
    
    if (is.null(lambda) | is.null(duration) | is.null(n_assessments))
      stop("input values are required for each of duration, n_assessments, and lambda")

    alternative <-
      match.arg(alternative, c('two.sided', 'one.sided')) #defaults to 'two.sided'
    
    if (!is.null(n)){
      n_2 <- n / lambda
      N <- n_2 + n
    }
    
    if (is.null(n))
      N <- NULL

    C <- (sig2.w*12*(n_assessments-1)/(n_assessments*(n_assessments+1)*duration^2)+sig2.b)

    n.body <- quote({
      # sp=(1+ratio)^2*(qnorm(1-alpha)+qnorm(powerinput))^2*C/(bet^2*ratio)
      return((1 + lambda) ^ 2 * (qnorm(
        1 -
          ifelse(alternative == "two.sided", sig.level / 2, sig.level)
      ) + qnorm(power)) ^ 2  * C/
          (delta ^ 2 * lambda))
    })
    
    power.body <- quote({
      d <- delta/sqrt(C*(1/n+1/n_2))
      if(alternative == "two.sided"){
        # pr=1-(pnorm(qnorm(1-alpha/2)-d)-pnorm(qnorm(alpha/2)-d))
        return(1-(pnorm(qnorm(1-sig.level / 2)-d)-pnorm(qnorm(sig.level/2)-d)))
      }
      if(alternative == "one.sided"){
        # pr=1-pnorm(qnorm(1-alpha)-d)
        return(1-pnorm(qnorm(1-sig.level)-d))
      }
    })
    
    if (is.null(N))
      N <- eval(n.body)
    else
      if (is.null(sig.level))
        sig.level <- uniroot(function(sig.level)
          eval(n.body) - N,
          c(1e-10, 1 - 1e-10),
          tol = tol,
          extendInt = "yes")$root
    else
      if (is.null(power))
        power <- eval(power.body)
    else
      if (is.null(delta))
        delta <- uniroot(
          function(delta)
            eval(n.body) - N,
          sqrt(sig2.e) * c(1e-7, 1e+7),
          tol = tol,
          extendInt = "downX"
        )$root
    else
      # Shouldn't happen
      stop("internal error", domain = NA)
    
    if (is.null(n))
      n <- lambda * N / (1 + lambda)
    
    n_2 <- N - n
    
    structure(
      list(
        N = n + n_2,
        n = c(n, n_2),
        delta = delta,
        duration = duration,
        n_assessments = n_assessments,
        sig2.b  = sig2.b  ,
        sig2.w   = sig2.w,
        sig.level = sig.level,
        power = power,
        alternative = alternative,
        note = "N is *total* sample size and n is sample size in *each* group",
        method = "Two-stage mixed effects model. Fitzmaurice, Laird, and Ware (2011)."
      ),
      class = "power.longtest"
    )
  }