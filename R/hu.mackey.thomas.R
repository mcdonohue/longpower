hu.mackey.thomas.power <- function(
  bet,
  sig0,
  rho,
  sig1,
  sig01,
  sig, 
  n_t,
  r
){
  if (input$time_length==1){
    C=0
    h=input$h
    for (i in 1:(n_t-2))
    {
      ti<-seq(0,i*h,h)
      ti_mean<-mean(ti)
      ti2_mean<-mean(ti^2)
      Ci<- exp(-i*h*r)*(1-exp(-h*r))*(i+1)*(sig^2*ti2_mean+(i+1)*sig0^2*(ti2_mean-ti_mean^2))/(sig^4+(i+1)*sig^2*(sig1^2*ti2_mean+sig0^2+2*sig01*ti_mean)+(i+1)^2*(sig0^2*sig1^2-sig01^2)*(ti2_mean-ti_mean^2))
      C<-C+Ci
    }
    
    t<-seq(0,(n_t-1)*h,h)
    t_mean<-mean(t)
    t2_mean<-mean(t^2)
    
    C<-C+exp(-(n_t-1)*h*r)*n_t*(sig^2*t2_mean+n_t*sig0^2*(t2_mean-t_mean^2))/(sig^4+n_t*sig^2*(sig1^2*t2_mean+sig0^2+2*sig01*t_mean)+n_t^2*(sig0^2*sig1^2-sig01^2)*(t2_mean-t_mean^2))
  }
  
  if (input$time_length==2){
    C=0
    length_char=input$gap_length
    t=as.numeric(strsplit(length_char,split=",")[[1]])
    n_t=length(t)
    # t_mean=mean(t)
    # t2_mean=mean(t^2)
    for (i in 1:(n_t-2)){
      ti=t[1:(i+1)]
      ti_mean<-mean(ti)
      ti2_mean<-mean(ti^2)
      Ci<-(exp(-r*t[i+1])-exp(-r*t[i+2]))*(i+1)*(sig^2*ti2_mean+(i+1)*sig0^2*(ti2_mean-ti_mean^2))/(sig^4+(i+1)*sig^2*(sig1^2*ti2_mean+sig0^2+2*sig01*ti_mean)+(i+1)^2*(sig0^2*sig1^2-sig01^2)*(ti2_mean-ti_mean^2))
      C<-C+Ci
    }
    t_mean<-mean(t)
    t2_mean<-mean(t^2)
    
    C<-C+exp(-r*t[n_t])*n_t*(sig^2*t2_mean+n_t*sig0^2*(t2_mean-t_mean^2))/(sig^4+n_t*sig^2*(sig1^2*t2_mean+sig0^2+2*sig01*t_mean)+n_t^2*(sig0^2*sig1^2-sig01^2)*(t2_mean-t_mean^2))
    
  }
  
  if (input$level==1)
  {
    alpha=input$oneside
    if (input$choice==2)
    {
      m1=ifelse(is.null(input$n_tre),100,input$n_tre)
      m2=ifelse(is.null(input$n_pla),100,input$n_pla)
      
      var_btre<-(m1+m2)/(C*m1*m2)
      d<-bet/sqrt(var_btre)
      pr=1-pnorm(qnorm(1-alpha)-d)
      pr=round(pr,4)
      # pr<-pwr.t.test(d=d,n=(m1+m2),sig.level=alpha,type="one.sample",alternative="two.sided")
      # return(paste("The RCRM statistical power is ", as.numeric(format(round(pr, 3), nsmall = 2)), sep = ""))
      return(c(pr,C,bet,m1,m2))
    }
    
    
    else{
      ratio=input$ratio
      powerinput=ifelse(is.null(input$power),0.8,input$power)
      
      sp=(1+ratio)^2*(qnorm(1-alpha)+qnorm(powerinput))^2/(bet^2*C*ratio)
      spsize=ceiling(sp)
      # return(paste("The total sample size needed is ", spsize, sep = ""))
      return(c(spsize,ratio,powerinput))
    }
  }
  if (input$level==2)
  {
    alpha=input$twoside
    if (input$choice==2)
    {
      m1=ifelse(is.null(input$n_tre),100,input$n_tre)
      m2=ifelse(is.null(input$n_pla),100,input$n_pla)
      
      var_btre<-(m1+m2)/(C*m1*m2)
      d<-bet/sqrt(var_btre)
      pr=1-(pnorm(qnorm(1-alpha/2)-d)-pnorm(qnorm(alpha/2)-d))
      pr=round(pr,4)
      # pr<-pwr.t.test(d=d,n=(m1+m2),sig.level=alpha,type="one.sample",alternative="two.sided")
      # return(paste("The RCRM statistical power is ", as.numeric(format(round(pr, 3), nsmall = 2)), sep = ""))
      return(c(pr,C,bet,m1,m2))
      
      
    }
    
    
    else{
      ratio=input$ratio
      powerinput=ifelse(is.null(input$power),0.8,input$power)
      
      sp=(1+ratio)^2*(qnorm(1-alpha/2)+qnorm(powerinput))^2/(bet^2*C*ratio)
      spsize=ceiling(sp)
      # return(paste("The total sample size needed is ", spsize, sep = ""))
      return(c(spsize,ratio,powerinput))
    }
  }
}