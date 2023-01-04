#' Random coefficient regression models (RCRM) sample size calculations
#'
#' This function computes sample size and power needed for the random coefficient 
#' regression models (RCRM) based on the formula from Hu, Mackey, and Thomas (2021). 
#' The RCRM assumes that the experimental and control arms have the same 
#' population baseline value. 
#' 
#' See Hu. Mackey, and Thomas (2021) for parameter details.
#'
#' @param N The total sample size. This formula can accommodate unbalanced
#' group allocation via \code{lambda}.
#' @param lambda allocation ratio (sample size group 1 divided by sample size group 2)
#' @param delta Effect size (absolute difference in rate of decline between tx and placebo)
#' @param t Vector of visit time points (including time 0)
#' @param sig2.i Variance of random intercept
#' @param cor.s.i Correlation between random intercept & slope
#' @param sig2.s Variance of random slope
#' @param sig2.e Variance of pure error
#' @param droprate Exponential dropout rate
#' @param sig.level type one error
#' @param p proportion vector for both groups; if i indexes visits, p[i] = the 
#' proportion whose last visit was at visit i (p sums to 1)
#' @param power power
#' @param alternative one- or two-sided test
#' @param tol	numerical tolerance used in root finding
#' @return One of the number of subject required per arm, the `power`, or
#' detectable effect size given `sig.level` and the other parameter estimates.
#'
#' @details See Equations (7) and (8) in Hu, Mackey, and Thomas (2021)
#' @author Monarch Shah
#' @seealso \code{\link{lmmpower}}, \code{\link{edland.linear.power}}, 
#' \code{\link{two.stage.me.power}}
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
t <- seq(0,1.5,0.25)
p <- c(rep(0, 6),1)
rcrm.power(delta=1.5, t=t, p=p, sig2.s = 24, sig2.e = 10, cor.s.i=0.5,
  sig.level=0.05, power = 0.80)
rcrm.power(N=360, t=t, sig2.s = 24, sig2.e = 10, cor.s.i=0.5,
  sig.level=0.05, power = 0.80)
rcrm.power(N=360, delta=1.5, t=t, sig2.s = 24, sig2.e = 10,
  cor.s.i=0.5, sig.level=0.05)

rcrm.power(delta=1.5, t=t, sig2.s = 24, sig2.e = 10, cor.s.i=0.5,
  sig.level=0.05, power = 0.80, alternative = 'one.sided')
rcrm.power(N=284, t=t, sig2.s = 24, sig2.e = 10, cor.s.i=0.5,
  sig.level=0.05, power = 0.80, alternative = 'one.sided')
rcrm.power(N=284, delta=1.5, t=t, sig2.s = 24, sig2.e = 10,
  cor.s.i=0.5, sig.level=0.05, alternative = 'one.sided')

rcrm.power <-
  function(N = NULL,
    delta = NULL,
    power = NULL,
    t = NULL,
    lambda = 1,
    sig2.i = 0,
    cor.s.i = NULL,
    sig2.s   = 0,
    sig2.e = NULL,
    sig.level = 0.05,
    p = NULL,
    alternative = c("two.sided", "one.sided"),
    tol = .Machine$double.eps^2)
  {
    # adapted from https://github.com/nan-hu-personal/RCRM_power_size/blob/master/app.R
    # https://rcrm-power-size.shinyapps.io/
    if (sum(sapply(list(N, delta,  power), is.null)) != 1)
      stop("exactly one of 'N', 'delta', and 'power' must be NULL")
    
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 >
        sig.level | sig.level > 1))
      stop("'sig.level' must be numeric in [0, 1]")
    
    if (!is.null(power) && !is.numeric(power) || any(0 >
        power | power > 1))
      stop("'power' must be numeric in [0, 1]")
    
    if (is.null(cor.s.i) |
        is.null(sig2.e) | is.null(t) | is.null(lambda) | is.null(p))
      stop("input values are required for each of t, p, lambda, cor.s.i and sig2.e")
    
    if (t[1] != 0)
      stop("t[1] must be 0")
    
    # Sum of p terms should be = 1
    if (round(sum(p)) != 1) {
      stop(cat("\n Sum of p terms not equal to 1, your sum of p is ", sum(p), "\n"))
    }
    
    # cor.s.i should be between -1 to 1
    if (cor.s.i < -1 | cor.s.i > 1) {
      stop(cat("\n cor.s.i should range from -1 to 1, you have specified: ", cor.s.i, "\n"))
    }
    
    # Negative delta to abs value
    if (delta < 0) {
      warning(cat(
        "\n Converting negative delta to absolute value, you have specified: ",
        delta,
        "\n"
      ))
    }
    
    # Variance terms should be > 0
    if (!any(sig2.i > 0, sig2.s > 0, sig2.e > 0)) {
      stop(cat("\n sig2.i, sig2.s or sig2.e should be > 0 \n"))
    }
    
    # Significance level should be between 0 to 1
    if (sig.level < 0 | sig.level > 1) {
      stop(cat(
        "\n Significance level should range from 0 to 1, you have specified: ",
        sig.level,
        "\n"
      ))
    }
    
    # Sum of p terms should be = 1
    if (round(sum(p)) != 1) {
      stop(cat("\n Sum of p terms not equal to 1, your sum of p is ", sum(p), "\n"))
    }
    
    # p values should be between 0 & 1
    if (any(sapply(p, function(x)
      x < 0 | x > 1))) {
      stop(cat("\n p terms should be between 0 to 1 \n"))
    }
    
    # Scan t and if baseline (0) is not mentioned give out warning message
    
    # t and p length should match
    if (length(t) != length(p)) {
      stop(cat("\n Length of t and p do not match \n"))
    }
    
    alternative <-
      match.arg(alternative, c('two.sided', 'one.sided')) #defaults to 'two.sided'
    
    if (!is.null(N)){
      # Based on lambda (Sample Size allocation ratio between experimental group and control group)
      # and n (total sample size), compute sample size in each group.
      n <- (N * lambda) / (lambda + 1)
      n_2 <- N / (lambda + 1)
    }
    
    cumsumt <- cumsum(t)
    cumsumtsq <- cumsum(t ** 2)
    k <- seq(0, length(t) - 1)
    
    tbar <- rep(NA, length(t))
    tbarsq <- rep(NA, length(t))
    
    num <- rep(NA, length(t))
    deno <- rep(NA, length(t))
    cstar <- rep(NA < length(t))
    
    for (i in 1:length(t)) {
      tbar[i] <- 1 / (k[i] + 1) * cumsumt[i]
      tbarsq[i] <- 1 / (k[i] + 1) * cumsumtsq[i]
      
      # Denominator
      deno[i] <- sig2.e ** 2 +
        (k[i] + 1) * sig2.e * (sig2.s * tbarsq[i] + sig2.i + 2 * cor.s.i * sqrt(sig2.i) * sqrt(sig2.s) * tbar[i]) +
        (k[i] + 1) ** 2 * sig2.i * sig2.s * (1 - cor.s.i ** 2) * (tbarsq[i] - tbar[i] * tbar[i])
      
      # Numerator
      num[i] <-
        (k[i] + 1) * (sig2.e * tbarsq[i] + (k[i] + 1) * sig2.i * (tbarsq[i] - tbar[i] * tbar[i]))
      
      cstar[i] <- p[i] * (num[i] / deno[i])
      
    }
        
    N.body <- quote({
      n_part1 <- (1 + lambda) ** 2 / (lambda * sum(cstar))
      n_part2 <- ((qnorm(1 - ifelse(alternative == "two.sided", sig.level / 2, sig.level)) + qnorm(power)) / abs(delta)) ** 2
      
      return(n_part1 * n_part2)
    })
    
    power.body <- quote({
      var_betahat <-
        (n + n_2) / (sum(cstar) * n * n_2)
      
      power_int <-
        abs(delta) / sqrt(var_betahat) - qnorm((1 - ifelse(alternative == "two.sided", sig.level / 2, sig.level)))
      
      return(pnorm(power_int))
    })
    
    if (is.null(N))
      N <- eval(N.body)
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
        p = p,
        sig2.i = sig2.i  ,
        sig2.s  = sig2.s  ,
        sig2.e   = sig2.e,
        cor.s.i = cor.s.i,
        sig.level = sig.level,
        power = power,
        alternative = alternative,
        note = "N is *total* sample size and n is sample size in *each* group",
        method = "Random Coefficients Regression Model. Hu, Mackey & Thomas (2021)"
      ),
      class = "power.longtest"
    )
  }