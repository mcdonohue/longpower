#' Title Random coefficient regression models (RCRM) sample size calculations
#' @description This function computes sample size and power needed for the random coefficient regression models (RCRM) based on the formula from Hu, Mackey, and Thomas (2021). The RCRM assumes that the experimental and control arms have the same population baseline value. See Hu. Mackey, and Thomas (2021) for parameter details.
#'
#' @param n total sample size
#' @param delta absolute value of slope difference between experimental and control arms
#' @param power power
#' @param t the observation times, including baseline 0
#' @param gamma sample size allocation ratio between the experimental group and control group (experimental over control) at baseline
#' @param rho correlation between random intercept and random slope
#' @param sig2.alpha random intercept variance
#' @param sig2.beta random slope variance
#' @param sig2 pure error variance
#' @param sig.level type one error
#' @param p proportion vector for both groups; if i indexes visits, p[i] = the proportion whose last visit was at visit i (p sums to 1)
#' @param alternative one- or two-sided test (must be one of "two.sided" or "one.sided")
#'
#' @details See Equations (7) and (8) in Hu, Mackey, and Thomas (2021)
#'
#' @author Monarch Shah
#'
#' @return One of the total number of subject required, or the `power`
#' @export
#'
#' @examples
#' # given n, compute power
#' example1 <-
#'  hu.mackey.thomas.linear.power(
#'   n = 368,
#'   delta = .33,
#'   t = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5),
#'   gamma = 1,
#'   rho = 0.07,
#'   sig2.alpha = 0.54,
#'   sig2.beta = 1.01,
#'   sig2 = 0.81,
#'   sig.level = .05, alternative = "one.sided",
#'   p = c(.0397, .0381, .0366, .0351, .2740, .2535, .2343, .0886)
#' )
#'
#' example1
#'
#' # given power, compute n
#' example2 <-
#'   hu.mackey.thomas.linear.power(
#'    power = .8,
#'    delta = .33,
#'    t = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5),
#'    gamma = 1,
#'    rho = 0.07,
#'    sig2.alpha = 0.54,
#'    sig2.beta = 1.01,
#'    sig2 = 0.81,
#'    sig.level = .05,
#'    p = c(.0397, .0381, .0366, .0351, .2740, .2535, .2343, .0886)
#'  )
#'
#'example2

hu.mackey.thomas.linear.power <-
  function(n = NULL,
           delta = NULL,
           power = NULL,
           t = NULL,
           gamma = 1,
           rho = NULL,
           sig2.alpha = NULL,
           sig2.beta = NULL,
           sig2 = NULL,
           sig.level = 0.05,
           p = NULL,
           alternative = "two.sided") {
    # Warnings & Error Message Handling
    # User should only input either n or power
    if (sum(is.null(n), is.null(power)) != 1) {
      stop(cat("\n Exactly one of 'n' or 'power' must be NULL \n"))
    }
    
    # Expected parameters are missing
    if (any(
      is.null(delta),
      is.null(t),
      is.null(gamma),
      is.null(rho),
      is.null(sig2.alpha),
      is.null(sig2.beta),
      is.null(sig2.beta),
      is.null(sig2),
      is.null(p)
    )) {
      stop(cat("\n One or more expected parameters are missing \n"))
    }
    
    # Power should be between 0 to 1
    if (!is.null(power)) {
      if (power < 0 | power > 1) {
        stop(cat("\n Power should range from 0 to 1, you have specified: ", power, "\n"))
      }
    }
    
    # Rho should be between -1 to 1
    if (rho < -1 | rho > 1) {
      stop(cat("\n Rho should range from -1 to 1, you have specified: ", rho, "\n"))
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
    if (!any(sig2.alpha > 0, sig2.beta > 0, sig2 > 0)) {
      stop(cat("\n sig2.alpha, sig2.beta or sig2 should be > 0 \n"))
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
    
    # Based on gamma (Sample Size allocation ratio between experimental group and control group)
    # and n (total sample size), compute sample size in each group.
    n_grp1 <- (n * gamma) / (gamma + 1)
    n_grp2 <- n / (gamma + 1)
    
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
      deno[i] <- sig2 ** 2 +
        (k[i] + 1) * sig2 * (sig2.beta * tbarsq[i] + sig2.alpha + 2 * rho * sqrt(sig2.alpha) * sqrt(sig2.beta) * tbar[i]) +
        (k[i] + 1) ** 2 * sig2.alpha * sig2.beta * (1 - rho ** 2) * (tbarsq[i] - tbar[i] * tbar[i])
      
      # Numerator
      num[i] <-
        (k[i] + 1) * (sig2 * tbarsq[i] + (k[i] + 1) * sig2.alpha * (tbarsq[i] - tbar[i] * tbar[i]))
      
      cstar[i] <- p[i] * (num[i] / deno[i])
      
    }
    
    # Alpha set up
    if (alternative == "two.sided") {
      alpha <- sig.level / 2
    }
    else if (alternative == "one.sided") {
      alpha <- sig.level
    }
    
    # Power calculations
    if (is.null(power)) {
      var_betahat <-
        (n_grp1 + n_grp2) / (sum(cstar) * n_grp1 * n_grp2)
      
      power_int <-
        abs(delta) / sqrt(var_betahat) - qnorm((1 - alpha))
      
      power <- pnorm(power_int)
      
      # return(power)
      
    }
    
    if (is.null(n)) {
      n_part1 <- (1 + gamma) ** 2 / (gamma * sum(cstar))
      n_part2 <- ((qnorm(1 - alpha) + qnorm(power)) / abs(delta)) ** 2
      n <- ceiling(n_part1 * n_part2)
      
      # return(n)
      
    }
    
    METHOD <- "Hu, Mackey, and Thomas (2021)"
    structure(
      list(
        n = n,
        delta = delta,
        power = power,
        t = t,
        gamma = gamma,
        rho = rho,
        sig2.alpha = sig2.alpha,
        sig2.beta = sig2.beta,
        sig2 = sig2,
        sig.level = sig.level,
        p = p,
        alternative = alternative,
        note = "n is *total* sample size",
        method = METHOD
      ),
      class = "power.longtest"
    )
    
  }
