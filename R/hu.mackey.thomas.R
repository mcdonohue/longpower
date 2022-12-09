#' Power for Random Coefficient Regression Model (RCRM)
#'
#' The RCRM theoretical results are summarized in a manuscript (Hu, Mackey, and
#' Thomas; in progress). The two-stage mixed effects model details can be found
#' in Section 8.4 of Fitzmaurice, Laird, and Ware (2011).
#'
#' This program was co-developed by Nan Hu (Genentech Inc., Biostatistics)
#' and a summer intern Zhe Qu (Ph.D. candidate at Tulane University). If you
#' have any technical questions or any ideas for improving functionality,
#' please contact Nan Hu (hu.nan@@gene.com).
#'
#' Hu, N., Mackey, H., & Thomas, R. (2021). Power and sample size for random coefficient regression models in randomized experiments with monotone missing data. Biometrical Journal, 63(4), 806-824.
#'
#' @param n sample size, group 1
#' @param lambda allocation ratio (= (sample size group 1)/(sample size group 2)
#' @param delta Effect size (absolute difference in rate of decline between tx and placebo)
#' @param t Vector of visit time points (including time 0)
#' @param sig2.int Variance of random intercept
#' @param rho Correlation between random intercept & slope
#' @param sig2.s Variance of random slope
#' @param sig2.e Variance of pure error
#' @param droprate Exponential dropout rate
#' @param sig.level type one error
#' @param power power
#' @param alternative one- or two-sided test
#' @param tol	not used (no root finding used in this implementation).
#' @return One of the number of subject required per arm, the `power`, or
#' detectable effect size given `sig.level` and the other parameter estimates.
#'
#' @author Nan Hu, Michael C. Donohue
#' @seealso \code{\link{lmmpower}}, \code{\link{edland.linear.power}}
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
#' hu.mackey.thomas.power(delta=1.5, t=t, sig2.s = 24, sig2.e = 10, rho=0.5, sig.level=0.05, power = 0.80)
#' hu.mackey.thomas.power(n=180, t=t, sig2.s = 24, sig2.e = 10, rho=0.5, sig.level=0.05, power = 0.80)
#' hu.mackey.thomas.power(n=180, delta=1.5, t=t, sig2.s = 24, sig2.e = 10, rho=0.5, sig.level=0.05)
#'
#' hu.mackey.thomas.power(delta=1.5, t=t, sig2.s = 24, sig2.e = 10, rho=0.5, sig.level=0.05, power = 0.80, alternative = 'one.sided')
#' hu.mackey.thomas.power(n=142, t=t, sig2.s = 24, sig2.e = 10, rho=0.5, sig.level=0.05, power = 0.80, alternative = 'one.sided')
#' hu.mackey.thomas.power(n=142, delta=1.5, t=t, sig2.s = 24, sig2.e = 10, rho=0.5, sig.level=0.05, alternative = 'one.sided')
#'
hu.mackey.thomas.power <-
  function(n = NULL,
    delta = NULL,
    power = NULL,
    t = NULL,
    lambda = 1,
    sig2.int = 0,
    rho = NULL,
    sig2.s   = 0,
    sig2.e = NULL,
    sig.level = 0.05,
    droprate = NULL,
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
    
    if (is.null(rho) |
        is.null(sig2.e) | is.null(t) | is.null(lambda))
      stop("input values are required for each of t, lambda, rho and sig2.e")
    
    if (t[1] != 0)
      stop("t[1] must be 0")
    
    if (is.null(droprate))
      droprate <- 0
    
    alternative <-
      match.arg(alternative, c('two.sided', 'one.sided')) #defaults to 'two.sided'
    
    n_t <- length(t)
    sig01 <- rho * sqrt(sig2.int) * sqrt(sig2.s)
    if (!is.null(n)){
      n_2 <- n / lambda
      N <- n_2 + n
    }
    
    if (is.null(n))
      N <- NULL
    
    C <- 0
    n_t <- length(t)
    # t_mean<-mean(t)
    # t2_mean<-mean(t^2)
    for (i in 1:(n_t - 2)) {
      ti <- t[1:(i + 1)]
      ti_mean <- mean(ti)
      ti2_mean <- mean(ti ^ 2)
      Ci <-
        (exp(-droprate * t[i + 1]) - exp(-droprate * t[i + 2])) * (i + 1) * (sig2.e * ti2_mean +
            (i + 1) * sig2.int * (ti2_mean - ti_mean ^ 2)) / (
              sig2.e ^ 2 + (i + 1) * sig2.e * (sig2.s * ti2_mean + sig2.int + 2 * sig01 *
                  ti_mean) + (i + 1) ^ 2 * (sig2.int * sig2.s - sig01 ^ 2) * (ti2_mean -
                      ti_mean ^ 2)
            )
      C <- C + Ci
    }
    t_mean <- mean(t)
    t2_mean <- mean(t ^ 2)
    
    C <-
      C + exp(-droprate * t[n_t]) * n_t * (sig2.e * t2_mean + n_t * sig2.int * (t2_mean -
          t_mean ^ 2)) / (
            sig2.e ^ 2 + n_t * sig2.e * (sig2.s * t2_mean + sig2.int + 2 * sig01 * t_mean) +
              n_t ^ 2 * (sig2.int * sig2.s - sig01 ^ 2) * (t2_mean - t_mean ^ 2)
          )
    
    n.body <- quote({
      return((1 + lambda) ^ 2 * (qnorm(
        1 -
          ifelse(alternative == "two.sided", sig.level / 2, sig.level)
      ) + qnorm(power)) ^ 2 /
          (delta ^ 2 * C * lambda))
    })
    
    power.body <- quote({
      var_btre <- (n + n_2) / (C * n * n_2)
      d <- delta / sqrt(var_btre)
      return(1 - (pnorm(qnorm(
        1 - ifelse(alternative == "two.sided", sig.level / 2, sig.level)
      ) - d) - pnorm(qnorm(sig.level / 2) -
          d)))
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
        t = t,
        droprate = droprate,
        sig2.int = sig2.int  ,
        sig2.s  = sig2.s  ,
        sig2.e   = sig2.e,
        rho = rho,
        sig.level = sig.level,
        power = power,
        alternative = alternative,
        note = "N is *total* sample size and n is sample size in *each* group",
        method = "Hu, Mackey & Thomas (2021). Biometrical Journal."
      ),
      class = "power.longtest"
    )
  }