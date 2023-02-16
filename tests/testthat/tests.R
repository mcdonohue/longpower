library(longpower)
fm1 <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm2 <- nlme::lme(Reaction ~ Days, random=~Days|Subject, sleepstudy)
fm3 <- nlme::lme(Reaction ~ Days, random=~1|Subject, sleepstudy)
trash <- capture.output(capture.output(fm4 <- gee::gee(Reaction ~ Days, id = Subject,
  data = sleepstudy,
  corstr = "exchangeable"), type = 'message'), type = 'output')

context("Power/sample size calculations")

test_that("Reproduce table of n per group for exchangeable correlation on page 29 of Diggle et al. (1994)", {
  n = 3
  t = c(0,2,5)
  u = list(u1 = t, u2 = rep(0,n))
  v = list(v1 = cbind(1,1,rep(0,n)),
           v2 = cbind(1,0,t))         
  rho = c(0.2, 0.5, 0.8)
  sigma2 = c(100, 200, 300)
  tab = outer(rho, sigma2, 
        Vectorize(function(rho, sigma2){
          ceiling(diggle.linear.power(
            d=0.5,
            t=t,
            sigma2=sigma2,
            R=rho,
            alternative="one.sided",
            power=0.80)$n[1])}))
  expect_identical(tab, 
    matrix(c(313, 625, 938,
             196, 391, 586,
             79,  157, 235), ncol = 3, byrow = TRUE)
  )
})

test_that("Reproduce Table 1 from Lu, Luo, & Chen (2008)", {
  phi1 <- c(rep(1, 6), 2, 2)
  phi2 <- c(1, 1, rep(2, 6))
  lambda <- c(1, 2, sqrt(1/2), 1/2, 1, 2, 1, 2)
  ztest <- ttest1 <- c()
  for(i in 1:8){
    Na <- (phi1[i] + lambda[i] * phi2[i])*(qnorm(0.05/2) + qnorm(1-0.90))^2*(0.5^-2)
    Nb <- Na/lambda[i]
    ztest <- c(ztest, Na + Nb)
    v <- Na + Nb - 2
    Na <- (phi1[i] + lambda[i] * phi2[i])*(qt(0.05/2, df = v) + qt(1-0.90, df = v))^2*(0.5^-2)
    Nb <- Na/lambda[i]
    ttest1 <- c(ttest1, Na + Nb)
  }

  tab <- do.call(cbind, list(phi1, phi2, lambda, ztest, ttest1))
  tab2 <- do.call(cbind, list(c(1, 1, 1, 1, 1, 1, 2, 2), 
    c(1, 1, 2, 2, 2, 2, 2, 2), 
    c(1, 2, 0.707106781186548, 0.5, 1, 2, 1, 2), 
    c(168.11876898305, 189.133615105931, 244.966998329937, 252.178153474575, 252.178153474575, 315.222691843219, 
      336.2375379661, 378.267230211862), c(170.147602556564, 191.157604171059, 246.982211011525, 254.192519905287, 
        254.192519905287, 317.23132678479, 338.244744929957, 380.272061544442)))
  expect_equal(tab, tab2, tolerance = 1e-03)
})

test_that("Reproduce Table 2 from Lu, Luo, & Chen (2008)", {
  tab <- c()
  for(J in c(2,4))
  for(aJ in (1:4)/10)
  for(p1J in c(0, c(1, 3, 5, 7, 9)/10)){
    rJ <- 1-aJ
    r <- seq(1, rJ, length = J)
    # p1J = p^(J-1)
    tab <- c(tab, power.mmrm.ar1(rho = p1J^(1/(J-1)), ra = r, sigmaa = 1,
      lambda = 1, times = 1:J,
      delta = 1, sig.level = 0.05, power = 0.80)$phi1)
  }
  expect_equal(matrix(tab, ncol = 6, byrow = TRUE), 
    matrix(c(1.11111111111111, 1.11, 1.10111111111111, 1.08333333333333, 1.05666666666667, 1.02111111111111, 
             1.25, 1.2475, 1.2275, 1.1875, 1.1275, 1.0475, 
             1.42857142857143, 1.42428571428571, 1.39, 1.32142857142857, 1.21857142857143, 1.08142857142857, 
             1.66666666666667, 1.66, 1.60666666666667, 1.5, 1.34, 1.12666666666667, 
             1.11111111111111, 1.10050206679506, 1.08280448732802, 1.06283003239702, 1.03996686151349, 1.01408550239776, 
             1.25, 1.2247445108099, 1.18392908471514, 1.13886229048648, 1.08796962333115, 1.03089517940421, 
             1.42857142857143, 1.38254161826624, 1.31065330422472, 1.23318260782795, 1.14702199890211, 1.0514263712482, 
             1.66666666666667, 1.59010470647863, 1.47488822527398, 1.35405269950478, 1.221989488974, 1.07728040119987),
           ncol = 6, byrow = TRUE), tolerance = 1e-03)
})

test_that("lmmpower (diggle)", {
  expect_equal(lmmpower(delta=1.5, t = seq(0,1.5,0.25),
  	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)$n[1], 207.310093300903, 
    tolerance = 1e-03)
  expect_equal(lmmpower(n=208, t = seq(0,1.5,0.25),
  	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)$delta, 1.49751028943272, 
    tolerance = 1e-03)
  expect_equal(lmmpower(beta = 5, pct.change = 0.30, t = seq(0,1.5,0.25),
  	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)$n[1], 207.310093300903, 
    tolerance = 1e-03)
  expect_equal(lmmpower(fm1, pct.change = 0.30, t = seq(0,9,1), power = 0.80)$n[1], 68.4699286748582, 
    tolerance = 1e-03)
  expect_equal(lmmpower(fm2, pct.change = 0.30, t = seq(0,9,1), power = 0.80)$n[1], 68.4693809183413, 
    tolerance = 1e-03)
  # random intercept only
  expect_equal(lmmpower(fm3, pct.change = 0.30, t = seq(0,9,1), power = 0.80)$n[1], 18.5332152601122, 
    tolerance = 1e-03)
  expect_equal(lmmpower(fm4, pct.change = 0.30, t = seq(0,9,1), power = 0.80)$n[1], 18.845000035132, 
    tolerance = 1e-03)
})

test_that("power.mmrm.ar1", {
  Orthodont$t.index <- as.numeric(factor(Orthodont$age, levels = c(8, 10, 12, 14)))

  fmOrth.corAR1 <- nlme::gls( distance ~ Sex * I(age - 11), 
    Orthodont,
    correlation = corAR1(form = ~ t.index | Subject),
    weights = varIdent(form = ~ 1 | age) )

  C <- corMatrix(fmOrth.corAR1$modelStruct$corStruct)[[1]]
  sigmaa <- fmOrth.corAR1$sigma *
            coef(fmOrth.corAR1$modelStruct$varStruct, unconstrained = FALSE)['14']
  ra <- seq(1,0.80,length=nrow(C))
  
  expect_equal(
    power.mmrm(N=100, Ra = C, ra = ra, sigmaa = sigmaa, power = 0.80)$delta,
    power.mmrm.ar1(N=100, rho = C[1,2], ra = ra, sigmaa = sigmaa, power = 0.80)$delta, 
    tolerance = 1e-03
  )
})

test_that("Reproduce Diggle et al page 29 using liu.liang.linear.power", {
  t = c(0,2,5)
  n = length(t)
  u = list(u1 = t, u2 = rep(0,n))
  v = list(v1 = cbind(1,1,t),
           v2 = cbind(1,0,t))
  rho = c(0.2, 0.5, 0.8)
  sigma2 = c(100, 200, 300)
  tab.ll = outer(rho, sigma2, 
       Vectorize(function(rho, sigma2){
         ceiling(liu.liang.linear.power(
           delta=0.5, u=u, v=v,
           sigma2=sigma2,
           R=rho, alternative="one.sided",
           power=0.80)$N/2)}))
  expect_equal(
    as.numeric(tab.ll),
    c(313, 196, 79, 625, 391, 157, 938, 586, 235), 
    tolerance = 1e-03
  )
})

test_that("Reproduce Diggle et al page 30 using liu.liang.linear.power", {
  n = 3
  u = list(u1 = rep(1,n), u2 = rep(0,n))
  v = list(v1 = rep(1,n),
           v2 = rep(1,n))
  rho = c(0.2, 0.5, 0.8)
  delta = c(20, 30, 40, 50)/100
  tab = outer(rho, delta, 
       Vectorize(function(rho, delta){
         ceiling(liu.liang.linear.power(
           delta=delta, u=u, v=v,
           sigma2=1,
           R=rho, alternative="one.sided",
           power=0.80)$N/2)}))
  expect_equal(
    as.numeric(tab),
    c(145, 207, 268, 65, 92, 120, 37, 52, 67, 24, 33, 43), 
    tolerance = 1e-03
  )
})

test_that("lmmpower (liuliang)", {
  meth <- 'liuliang'
  expect_equal(lmmpower(delta=1.5, t = seq(0,1.5,0.25),
    sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), 
    power = 0.80, method = meth)$n[1], 207.310093300903, 
    tolerance = 1e-03)
  expect_equal(lmmpower(n=208, t = seq(0,1.5,0.25),
    sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), 
    power = 0.80, method = meth)$delta, 1.49751028943272, 
    tolerance = 1e-03)
  expect_equal(lmmpower(beta = 5, pct.change = 0.30, t = seq(0,1.5,0.25),
    sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), 
    power = 0.80, method = meth)$n[1], 207.310093300903, 
    tolerance = 1e-03)
  expect_equal(lmmpower(fm1, pct.change = 0.30, t = seq(0,9,1), 
    power = 0.80, method = meth)$n[1], 68.4699286748582, 
    tolerance = 1e-03)
  expect_equal(lmmpower(fm2, pct.change = 0.30, t = seq(0,9,1), 
    power = 0.80, method = meth)$n[1], 68.4693809183413, 
    tolerance = 1e-03)
  # random intercept only
  expect_equal(lmmpower(fm3, pct.change = 0.30, t = seq(0,9,1), 
    power = 0.80, method = meth)$n[1], 18.5332152601122, 
    tolerance = 1e-03)
  expect_equal(lmmpower(fm4, pct.change = 0.30, t = seq(0,9,1), 
    power = 0.80, method = meth)$n[1], 18.845000035132, 
    tolerance = 1e-03)
})

test_that("lmmpower (edland)", {
  meth <- 'edland'
  expect_equal(lmmpower(delta=1.5, t = seq(0,1.5,0.25),
    sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), 
    power = 0.80, method = meth)$n[1], 207.310093300903, 
    tolerance = 1e-03)
  expect_equal(lmmpower(n=208, t = seq(0,1.5,0.25),
    sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), 
    power = 0.80, method = meth)$delta, 1.49751028943272, 
    tolerance = 1e-03)
  expect_equal(lmmpower(beta = 5, pct.change = 0.30, t = seq(0,1.5,0.25),
    sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), 
    power = 0.80, method = meth)$n[1], 207.310093300903, 
    tolerance = 1e-03)
  expect_equal(lmmpower(fm1, pct.change = 0.30, t = seq(0,9,1), 
    power = 0.80, method = meth)$n[1], 68.4699286748582, 
    tolerance = 1e-03)
  expect_equal(lmmpower(fm2, pct.change = 0.30, t = seq(0,9,1), 
    power = 0.80, method = meth)$n[1], 68.4693809183413, 
    tolerance = 1e-03)
  # random intercept only
  expect_equal(lmmpower(fm3, pct.change = 0.30, t = seq(0,9,1), 
    power = 0.80, method = meth)$n[1], 18.5332152601122, 
    tolerance = 1e-03)
})

test_that("cprm.power", {
  meth <- 'cprm'
  # Table 3 in Zhao and Edland (2023)
  expect_equal(cprm.power(n = 80,t = 0:6*.25,delta = 3.488879,sig2.int = 55.45902,
                          sig.b0b1 = 13.58123,sig2.s = 15.44793,
                          sig2.e = 13.63768)$power, 0.7999999, 
               tolerance = 1e-03)
  expect_equal(cprm.power(power = 0.8,t = 0:6*.25,delta = 3.488879,sig2.int = 55.45902,
                          sig.b0b1 = 13.58123,sig2.s = 15.44793,
                          sig2.e = 13.63768)$n[1], 80.00001,
               tolerance = 1e-03)
  expect_equal(cprm.power(n = 80, power = 0.8,t = 0:6*.25,sig2.int = 55.45902,
                          sig.b0b1 = 13.58123,sig2.s = 15.44793,
                          sig2.e = 13.63768)$delta, 3.488879,
               tolerance = 1e-03)
  # Figure 3 in Zhao and Edland (2023)
  expect_equal(cprm.power(n = 200, delta = 0.25*4.0579*6*0.25,t = 0:6*.25,sig2.int = 55.45902,
                          sig.b0b1 = 13.58123,sig2.s = 15.44793,
                          sig2.e = 13.63768)$power, 0.4888692,
               tolerance = 1e-03)
  expect_equal(cprm.power(n = 400, delta = 0.25*4.0579*6*0.25,t = 0:6*.25,sig2.int = 55.45902,
                          sig.b0b1 = 13.58123,sig2.s = 15.44793,
                          sig2.e = 13.63768)$power, 0.7800555,
               tolerance = 1e-03)
  # random intercept only
  expect_equal(cprm.power(n = 80,t = 0:6*.25,delta = 1.5,sig2.int = 20,
                          sig.b0b1 = 0,sig2.s = 0,
                          sig2.e = 10)$power, 0.5640936, 
               tolerance = 1e-03)
})

test_that("lmmpower (hu)", {
  meth <- 'hu'
  expect_equal(lmmpower(delta=1.5, t = seq(0,1.5,0.25),
    sig2.i = 55, sig2.s = 24, sig2.e = 10, cor.s.i=0.8, p=c(rep(0,6),1), 
    power = 0.80, method = meth)$n[1], 135.4827, 
    tolerance = 1e-03)
  expect_equal(lmmpower(n=208, t = seq(0,1.5,0.25),
    sig2.i = 55, sig2.s = 24, sig2.e = 10, cor.s.i=0.8, p=c(rep(0,6),1), 
    power = 0.80, method = meth)$delta, 1.210603, 
    tolerance = 1e-03)
  expect_equal(lmmpower(beta = 5, pct.change = 0.30, t = seq(0,1.5,0.25),
    sig2.i = 55, sig2.s = 24, sig2.e = 10, cor.s.i=0.8, p=c(rep(0,6),1), 
    power = 0.80, method = meth)$n[1], 135.4827, 
    tolerance = 1e-03)
  expect_equal(lmmpower(fm1, pct.change = 0.30, t = seq(0,9,1), p=c(rep(0,9),1), 
    power = 0.80, method = meth)$n[1], 67.1745, 
    tolerance = 1e-03)
  expect_equal(lmmpower(fm2, pct.change = 0.30, t = seq(0,9,1), p=c(rep(0,9),1), 
    power = 0.80, method = meth)$n[1], 67.17401, 
    tolerance = 1e-03)
  # random intercept only
  expect_equal(lmmpower(fm3, pct.change = 0.30, t = seq(0,9,1), p=c(rep(0,9),1),  
    power = 0.80, method = meth)$n[1], 15.9781, 
    tolerance = 1e-03)
})

test_that("power.mmrm.ar1", {
  Orthodont$t.index <- as.numeric(factor(Orthodont$age, levels = c(8, 10, 12, 14)))
  
  fmOrth.corAR1 <- nlme::gls( distance ~ Sex * I(age - 11), 
    Orthodont,
    correlation = corAR1(form = ~ t.index | Subject),
    weights = varIdent(form = ~ 1 | age) )
  
  C <- corMatrix(fmOrth.corAR1$modelStruct$corStruct)[[1]]
  sigmaa <- fmOrth.corAR1$sigma *
    coef(fmOrth.corAR1$modelStruct$varStruct, unconstrained = FALSE)['14']
  ra <- seq(1,0.80,length=nrow(C))
  
  expect_equal(
    power.mmrm(N=100, Ra = C, ra = ra, sigmaa = sigmaa, power = 0.80)$delta,
    power.mmrm.ar1(N=100, rho = C[1,2], ra = ra, sigmaa = sigmaa, power = 0.80)$delta, 
    tolerance = 1e-03
  )
})

test_that("hu.mackey.thomas.linear.power vs original", {
  t <- seq(0,1.5,0.25)
  p <- t/sum(t)
  expect_equal(
    hu.mackey.thomas.linear.power(delta=1.5, t=t, sig2.s = 24, sig2.e = 10, cor.s.i=0.5, p=p, power = 0.80)$N,
    430.3265, 
    tolerance = 1e-03
  )
  expect_equal(
    hu.mackey.thomas.linear.power(n=216, t=t, sig2.s = 24, sig2.e = 10, cor.s.i=0.5, p=p, power = 0.80)$delta,
    1.497092, 
    tolerance = 1e-03
  )
  expect_equal(
    hu.mackey.thomas.linear.power(n=216, delta=1.5, t=t, sig2.s = 24, sig2.e = 10, cor.s.i=0.5, p=p)$power,
    0.8015201, 
    tolerance = 1e-03
  )
  expect_equal(
    hu.mackey.thomas.linear.power(delta=1.5, t=t, sig2.s = 24, sig2.e = 10, cor.s.i=0.5, p=p, power = 0.80, alternative = 'one.sided')$N,
    338.9679, 
    tolerance = 1e-04
  )
  expect_equal(
    hu.mackey.thomas.linear.power(n=170, t=t, sig2.s = 24, sig2.e = 10, cor.s.i=0.5, p=p, power = 0.80, alternative = 'one.sided')$delta,
    1.497722, 
    tolerance = 1e-03
  )
  expect_equal(
    hu.mackey.thomas.linear.power(n=170, delta=1.5, t=t, sig2.s = 24, sig2.e = 10, cor.s.i=0.5, p=p, alternative = 'one.sided')$power,
    0.8010573, 
    tolerance = 1e-03
  )
  
  #' Random coefficient regression models (RCRM) sample size calculations
  #' @param n total sample size (NOTE: different from longpower convention)
  #' @param gamma sample size allocation ratio between the experimental group and control group (experimental over control) at baseline
  #' @param rho correlation between random intercept and random slope
  #' @param sig2.alpha random intercept variance
  #' @param sig2.beta random slope variance
  #' @param sig2 pure error variance

  hu.mackey.thomas.linear.power.original <-
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
  
  example1_original <-
   hu.mackey.thomas.linear.power.original(
    n = 368,
    delta = .33,
    t = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5),
    gamma = 1,
    rho = 0.07,
    sig2.alpha = 0.54,
    sig2.beta = 1.01,
    sig2 = 0.81,
    p = c(.0397, .0381, .0366, .0351, .2740, .2535, .2343, .0886),
    alternative = "one.sided"
  )

  example1 <- hu.mackey.thomas.linear.power(n=368/2, delta = 0.33, 
    t=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5), cor.s.i = 0.07,
    sig2.i = 0.54, sig2.s = 1.01, sig2.e = 0.81, 
    p = c(.0397, .0381, .0366, .0351, .2740, .2535, .2343, .0886), 
    alternative = "one.sided")
  
  expect_equal(example1_original$power, example1$power)
  
  #' given power, compute n
  example2_orignal <- hu.mackey.thomas.linear.power.original(
     power = .8,
     delta = .33,
     t = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5),
     gamma = 1,
     rho = 0.07,
     sig2.alpha = 0.54,
     sig2.beta = 1.01,
     sig2 = 0.81,
     sig.level = .05,
     p = c(.0397, .0381, .0366, .0351, .2740, .2535, .2343, .0886)
   )
  
  example2 <- hu.mackey.thomas.linear.power(
    power = .8,
    delta = .33,
    t = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5),
    lambda = 1,
    cor.s.i = 0.07,
    sig2.i = 0.54,
    sig2.s = 1.01,
    sig2.e = 0.81,
    p = c(.0397, .0381, .0366, .0351, .2740, .2535, .2343, .0886)
  )
  
  expect_equal(example2_orignal$n, ceiling(example2$N))
})
