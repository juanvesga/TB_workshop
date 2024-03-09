## scale_up function
# Intervention scaling function
scale_up<- function (t, state, parameters,t.interv,parameters_old,fx) {
  
  scale <- min((t-t.interv[1])/(t.interv[2]-t.interv[1]),1); 
  if (scale<0) 
  {
    scale=0
  }
  
  pars_scaled <- parameters;
  
  pars_scaled <- parameters_old + scale*(parameters-parameters_old)
  
  return(fx(t, state,pars_scaled))
}


## int function

get_intervention<- function (sfin, params_new, params_old,times_new,
                             t.interv, fx_scale,fx_basic, int_name, data_stub) {
  
  # Starting conditions
  xstart <- c(U = sfin$U, 
              L = sfin$L,  
              I = sfin$I, 
              R = sfin$R,
              Ux = sfin$Ux, 
              Lx = sfin$Lx,
              Ix = sfin$Ix, 
              Rx = sfin$Rx,
              Up = sfin$Up, 
              Lp = sfin$Lp,  
              Ip = sfin$Ip, 
              Rp = sfin$Rp,
              Incidence= sfin$Incidence,
              Irecent=   sfin$Irecent,  
              Iremote=   sfin$Iremote,
              Incidence_x= sfin$Incidence_x,
              Irecent_x=   sfin$Irecent_x,  
              Iremote_x=   sfin$Iremote_x,
              Incidence_p = sfin$Incidence_p,
              Irecent_p=   sfin$Irecent_p,  
              Iremote_p=   sfin$Iremote_p) 

  #Select type of function
  if (sum(is.na(t.interv))>0)
  {
    fx<-fx_basic
  }  else {
    fx<-fx_scale
  }

  #Run the model
  out <- as.data.frame(ode(y = xstart, times = times_new, 
                           func = fx, parms = params_new))  #
                       
  # Model output
  N           <- out$U+out$L+out$I+out$R  
  prev        <- out$I/N
  rate.inc    <- 1e5*(diff(out$Incidence)/N[1:length(N)-1])
  fr.remo     <- diff(out$Iremote)/diff(out$Incidence)

  Np            <- out$Up+out$Lp+out$Ip+out$Rp  
  prev_p        <- out$Ip/Np
  rate.inc_p    <- 1e5*(diff(out$Incidence_p)/Np[1:length(Np)-1])
  fr.remo_p     <- diff(out$Iremote_p)/diff(out$Incidence_p)
  
  Nx            <- out$Ux+out$Lx+out$Ix+out$Rx  
  prev_x        <- out$Ix/Nx
  rate.inc_x    <- 1e5*(diff(out$Incidence_x)/Nx[1:length(Nx)-1])
  fr.remo_x     <- diff(out$Iremote_x)/diff(out$Incidence_x)

  time         <- out$time[1:length(out$time)-1]
  

  dat         <- data.frame(
    Years=time+(2024-400), 
    Incidence=rate.inc,
    Incidence_p=rate.inc_p,
    Incidence_x=rate.inc_x,
    prev=prev[c(2:length(N))],
    prev_x=prev_x[c(2:length(N))],
    prev_p=prev_p[c(2:length(N))],
    N=N[c(2:length(N))],
    Nx=Nx[c(2:length(Nx))],
    Np=Np[c(2:length(Np))])
  
  
  
  dat$Sim     <- int_name
  
  # If it is a first run, nothing to append 
  if (sum(is.na(data_stub))>0)
  {
    data<-dat
  }
  else # Append previous runs
  {
    data  <-rbind(data_stub, dat)
  }
  
  remote<-fr.remo  # rename a variable with the fraction of remote incidence 
  remote_p<-fr.remo_p  # rename a variable with the fraction of remote incidence 
  remote_x<-fr.remo_x  # rename a variable with the fraction of remote incidence 
  
  titl  <-int_name
  
  # Create our plot
  p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence, col=Sim))
  p1<-p + 
    geom_line(size=1.2) +
    ggtitle ('TB Incidence') +
    theme_bw() + ylab('Rate per 100,000 pop')+
    ylim(0,max(data$Incidence))

  p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence_p, col=Sim))
  p2<-p + 
    geom_line(size=1.2) +
    ggtitle ('TB Incidence in prisons') +
    theme_bw() + ylab('Rate per 100,000 imprisoned')+
    ylim(0,max(data$Incidence_p))

  p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence_x, col=Sim))
  p3<-p + 
    geom_line(size=1.2) +
    ggtitle ('TB Incidence among ex-prisoners') +
    theme_bw() + ylab('Rate per 100,000 ex-prisoners')+
    ylim(0,max(data$Incidence_x))

  p<- ggplot(data=data, mapping = aes(x=Years, y=Np*1e5, col=Sim))
  p4<-p + 
    geom_line(size=1.2) +
    ggtitle ('Imprisioned population') +
    theme_bw() + ylab('Population per 100k')+
    ylim(0,max(data$Np*1e5))
  
  p<- ggplot(data=data, mapping = aes(x=Years, y=Nx*1e5, col=Sim))
  p5<-p + 
    geom_line(size=1.2) +
    ggtitle ('Ex-prisoner population') +
    theme_bw() + ylab('Population per 100k')+
    ylim(0,max(data$Nx*1e5))

  df1<- data.frame(
    Years=time+(2024-400), 
    Incidence=rate.inc)
  df1$group<-c("Community")
  
  df2<- data.frame(
    Years=time+(2024-400), 
    Incidence=rate.inc_p)
  df2$group="Prisoners"
  
  df3<- data.frame(
    Years=time+(2024-400), 
    Incidence=rate.inc_x)
  df3$group="Ex-Prisoners"
  
  df<-rbind(df1,df2,df3)
  
  
  p<- ggplot(df, aes(x=Years, y=Incidence, color=group))
  p6<-p + 
    geom_line(size=1.2) +
    ggtitle ('TB Incidence') +
    theme_bw() + ylab('Rate per 100,000')+
    ylim(0,max(data$Incidence_p))
  
 
  df1<- data.frame(
    Years=time+(2024-400), 
    Prevalence=dat$prev*1e2)
  df1$group<-c("Community")
  
  df2<- data.frame(
    Years=time+(2024-400), 
    Prevalence=dat$prev_p*1e2)
  df2$group="Prisoners"
  
  df3<- data.frame(
    Years=time+(2024-400), 
    Prevalence=dat$prev_x*1e2)
  df3$group="Ex-Prisoners"
  
  df<-rbind(df1,df2,df3)
  
  
  p<- ggplot(df, aes(x=Years, y=Prevalence, color=group))
  p7<-p + 
    geom_line(size=1.2) +
    ggtitle ('TB Prevalence') +
    theme_bw() + ylab('%')+
    ylim(0,max(data$Incidence_p))
  
  

  # Pie chart of remote vs recent incidence 
  df <- data.frame(
    Source = c("Recent", "Remote"),
    value  = c(1-tail(remote,1),tail(remote,1))
  )
  df2 <- data.frame(
    Source = c("Recent", "Remote"),
    value  = c(1-tail(remote_p,1),tail(remote_p,1))
  )
  df3 <- data.frame(
    Source = c("Recent", "Remote"),
    value  = c(1-tail(remote_x,1),tail(remote_x,1))
  )
  
  mycols <- c("#0073C2FF", "#EFC000FF")
  pie1<- ggplot(df, aes(x="", y=value, fill=Source))+
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0)+
    scale_fill_manual(values = mycols) +
    theme_void()+
    ggtitle (titl) 

  pie2<- ggplot(df2, aes(x="", y=value, fill=Source))+
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0)+
    scale_fill_manual(values = mycols) +
    theme_void()+
    ggtitle (titl) 
  
  pie3<- ggplot(df3, aes(x="", y=value, fill=Source))+
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0)+
    scale_fill_manual(values = mycols) +
    theme_void()+
    ggtitle (titl) 
  
  output<-list("out"=out, 
               "inc_c"=p1,
               "inc_p"=p2,
               "inc_x"=p3,
               "pop_p"=p4,
               "pop_x"=p5, 
               "pie1"=pie1,
               "pie2"=pie2, 
               "pie3"=pie3, 
               "inc_all"=p6, 
               "prev_all"=p7, 
               "data"=data)
  
  # Spit-out results
  return(output)
  
}