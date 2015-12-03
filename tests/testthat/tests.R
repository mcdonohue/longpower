library(longpower)
context("Power/sample size calculations")

test_that("Reproduce table for exchangeable correlation on page 29 of Diggle et al. (1994)", {

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
            power=0.80)$n)}))
  expect_identical(tab, 
    matrix(c(313, 625, 938,
                196, 391, 586,
                79, 157, 235), ncol = 3, byrow = TRUE)
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
  expect_equal(tab, tab2)
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
           ncol = 6, byrow = TRUE))
})

test_that("lmmpower", {
  expect_equal(lmmpower(delta=1.5, t = seq(0,1.5,0.25),
  	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)$n, 207.310093300903)
  expect_equal(lmmpower(n=208, t = seq(0,1.5,0.25),
  	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)$delta, 1.49748557462127)
  expect_equal(lmmpower(beta = 5, pct.change = 0.30, t = seq(0,1.5,0.25),
  	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)$n, 207.310093300903)

  library(lme4)
  fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
  expect_equal(lmmpower(fm1, pct.change = 0.30, t = seq(0,9,1), power = 0.80)$n, 68.4699256666055)

  library(nlme)
  fm2 <- lme(Reaction ~ Days, random=~Days|Subject, sleepstudy)
  expect_equal(lmmpower(fm2, pct.change = 0.30, t = seq(0,9,1), power = 0.80)$n, 68.4693809183413)

  # random intercept only
  fm3 <- lme(Reaction ~ Days, random=~1|Subject, sleepstudy)
  expect_equal(lmmpower(fm3, pct.change = 0.30, t = seq(0,9,1), power = 0.80)$n, 18.5332152601122)

  library(gee)
  trash <- capture.output(capture.output(fm4 <- gee(Reaction ~ Days, id = Subject,
              data = sleepstudy,
              corstr = "exchangeable"), type = 'message'), type = 'output')
  expect_equal(lmmpower(fm4, pct.change = 0.30, t = seq(0,9,1), power = 0.80)$n, 18.845000035132)
})