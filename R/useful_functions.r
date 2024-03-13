## scale_up function
# Intervention scaling function
scale_up<- function (t, state, parameters,t.interv,parameters_old,fx) {
  
  scale <- min((t-t.interv[1])/(t.interv[2]-t.interv[1]),1); 
  if (scale<0) 
  {
    scale=0
  }
  
  pars_scaled <- parameters
  
  
  pars_scaled <- parameters_old + scale*(parameters-parameters_old)
  
  return(fx(t, state,pars_scaled))
}


## int function

get_intervention<- function (sfin, params_new, params_old,times_new,
                             t.interv, int_name, data_stub=NA) {
  
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
    fx<-TB_prison
  }  else {
    
    fx<-function(t, state, parameters) scale_up(t, 
                                                state, 
                                                parameters,
                                                t.interv,
                                                params_old,
                                                TB_prison)
    
  }
  
  #Run the model
  out <- as.data.frame(ode(y = xstart, 
                           times = times_new, 
                           func = fx, 
                           parms = params_new))  #
  
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
  
  levs<-c("Baseline",
          "Diagnosis",
          "Treatment",
          "Fast access",
          "Prevention",
          "Fast access X",
          "Prevention X")
  
  data$Sim<-factor(data$Sim, levels=levs)
  remote<-fr.remo  # rename a variable with the fraction of remote incidence 
  remote_p<-fr.remo_p  # rename a variable with the fraction of remote incidence 
  remote_x<-fr.remo_x  # rename a variable with the fraction of remote incidence 
  
  titl  <-int_name
  
  
  colors <- c("Community" = "gold", 
              "Ex-prisoners" = "lightseagreen",
              "Prisoners" = "grey43")
  
  colors_scen<-c("Baseline"="#c23d3d",
          "Diagnosis"="#5342f0",
          "Treatment"="#ff6600",
          "Fast access"="#055770",
          "Prevention"="#d3b016",
          "Fast access X" = "#3aaf1d",
          "Prevention X"="#e448e4")
  
  # Create our plot
  p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence))
  pb1<-p + 
    geom_line(linewidth=1.2, color="firebrick") +
    geom_point(aes(x=2023,y=250),shape = 15, color = "black", size = 3) +
    ggtitle ('TB Incidence (Community)') +
    theme_bw() + ylab('Rate per 100k')+
    ylim(c(0,1000))+
    xlim(c(1980,2024))
  
  p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence_p))
  pb2<-p + 
    geom_line(linewidth=1.2,, color="firebrick") +
    geom_point(aes(x=2023,y=990),shape = 15, color = "black", size = 3) +
    ggtitle ('TB Incidence (Prisons)') +
    theme_bw() + ylab('Rate per 100k')+
     ylim(c(0,1000))+
    xlim(c(1980,2024))
  
  p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence_x))
  pb3<-p + 
    geom_line(linewidth=1.2) +
    ggtitle ('TB Incidence (Ex-prisoners)') +
    theme_bw() + ylab('Rate per 100k')+
    ylim(0,max(data$Incidence_x))
  
  
  
  p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence, col=Sim))
  p1<-p + 
    geom_line(linewidth=1.2) +
    ggtitle ('TB Incidence (Community)') +
    theme_bw() + ylab('Rate per 100k')+
    theme(legend.title=element_blank(), legend.text=element_text(size=7))+ 
    ylim(0,max(data$Incidence))+
    xlim(c(2024,2035))+
    scale_color_manual(values = colors_scen) 
  
  p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence_p, col=Sim))
  p2<-p + 
    geom_line(linewidth=1.2) +
    ggtitle ('TB Incidence (Prisons)') +
    theme_bw() + ylab('Rate per 100k')+
    theme(legend.title=element_blank(), legend.text=element_text(size=7))+ 
    ylim(0,max(data$Incidence_p))+
    xlim(c(2024,2035))+
    scale_color_manual(values = colors_scen) 
  
  p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence_x, col=Sim))
  p3<-p + 
    geom_line(linewidth=1.2) +
    ggtitle ('TB Incidence (Ex-prisoners)') +
    theme_bw() + ylab('Rate per 100k')+
    theme(legend.title=element_blank(), legend.text=element_text(size=7))+ 
    ylim(0,max(data$Incidence_x))+
    xlim(c(2024,2035))+
    scale_color_manual(values = colors_scen) 
  
  p<- ggplot(data=data, mapping = aes(x=Years, y=Np*1e5, col=Sim))
  p4<-p + 
    geom_line(linewidth=1.2) +
    ggtitle ('Imprisioned population') +
    theme_bw() + ylab('Population per 100k')+
    ylim(0,max(data$Np*1e5))
  
  p<- ggplot(data=data, mapping = aes(x=Years, y=Nx*1e5, col=Sim))
  p5<-p + 
    geom_line(linewidth=1.2) +
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
  df3$group="Ex-prisoners"
  
  df<-rbind(df1,df2,df3)
  

  
  p<- ggplot(df, aes(x=Years, y=Incidence, color=group))
  p6<-p + 
    geom_line(linewidth=1.2) +
    ggtitle ('TB Incidence') +
    theme_bw() + ylab('Rate per 100,000')+
    scale_color_manual(values = colors) +
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
    geom_line(linewidth=1.2) +
    ggtitle ('TB Prevalence') +
    theme_bw() + ylab('%')+
    ylim(0,max(data$Incidence_p))
  
  p8=NA
  
  if(int_name%in%levs){
 

    df<-data[data$Sim==int_name,] 
    idend<-25
    id0<-1

    reductions<-data.frame(
      y=c(
        1e2*(1-(df$Incidence[idend]/df$Incidence[id0])),
        1e2*(1-(df$Incidence_x[idend]/df$Incidence_x[id0])),
        1e2*(1-(df$Incidence_p[idend]/df$Incidence_p[id0]))
      ),
      Group=c("Community", "Ex-prisoners","Prisoners"))
    
   
    p8<-ggplot(reductions, aes(x=Group, y=y, fill=Group))+
      geom_bar(stat = "identity")+
      ggtitle ('TB Incidence reductions in 10 years') +
      ylim(c(0,75))+
      theme_bw() + ylab('%')+xlab("")+
      scale_fill_manual(values = colors) +
      theme(
        legend.position = "none",
        panel.background = element_blank()
      )
  }
  
  
  output<-list("out"=out, 
               "incbase_c"=pb1,
               "incbase_p"=pb2,
               "incbase_x"=pb3,
               "inc_c"=p1,
               "inc_p"=p2,
               "inc_x"=p3,
               "pop_p"=p4,
               "pop_x"=p5, 
               "inc_all"=p6, 
               "prev_all"=p7,
               "reduc"=p8,
               "data"=data)
  
  # Spit-out results
  return(output)
  
}