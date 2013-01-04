liu.liang.linear.power <-
function(delta, u, v, sigma2=1, R, 
         sig.level=0.05, power=0.80, 
         Pi = rep(1/length(u),length(u)),
         alternative = c("two.sided", "one.sided"))
{
   alternative = match.arg(alternative)
   sig.level = ifelse(alternative=="two.sided", sig.level/2, sig.level)
   if(is.null(dim(R))) R = matrix(R, length(u[[1]]), length(u[[1]])) + diag(1-R,length(u[[1]]))

   Ipl = 0
   for(i in 1:length(u)) 
      Ipl = Ipl + Pi[i]*t(u[[i]])%*%solve(R)%*%v[[i]]
   Ipl = Ipl/sigma2   
      
   Ill = 0
   for(i in 1:length(u)) 
      Ill = Ill + Pi[i]*t(v[[i]])%*%solve(R)%*%v[[i]]
   Ill = Ill/sigma2

   Sigma1 = 0
   for(i in 1:length(u)) 
      Sigma1 = Sigma1 + Pi[i]*(t(u[[i]])-Ipl%*%solve(Ill)%*%t(v[[i]]))%*%solve(R)%*%
                              (u[[i]]-v[[i]]%*%solve(Ill)%*%t(Ipl))                           
   Sigma1 = Sigma1/sigma2

   n = ceiling((qnorm(1-sig.level) + qnorm(power))^2/(delta%*%Sigma1%*%delta)[1,1])
   NOTE <- "n is number in *each* group"
   METHOD <- "Longitudinal linear model power calculation (Liu & Liang, 1997)"
   structure(list(n = n, delta = delta, sigma2 = sigma2, R = R, 
        sig.level = sig.level, power = power, alternative = alternative, note = NOTE, 
        method = METHOD), class = "power.longtest")

}

