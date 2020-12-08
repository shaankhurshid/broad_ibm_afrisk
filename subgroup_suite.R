# Script to explore coefficients in subgroups

## Source helper functions
# A. KM function
kmdec=function(dec.num,dec.name, datain, adm.cens){
  stopped=0
  data.sub=datain[datain[,dec.name]==dec.num,]
  if (sum(data.sub$out)>1){
    avsurv=survfit(Surv(tvar,out) ~ 1, data=datain[datain[,dec.name]==dec.num,], error="g")
    avsurv.est=ifelse(min(avsurv$time)<=adm.cens,avsurv$surv[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],1)
    
    avsurv.stderr=ifelse(min(avsurv$time)<=adm.cens,avsurv$std.err[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],0)
    avsurv.stderr=avsurv.stderr*avsurv.est
    
    avsurv.num=ifelse(min(avsurv$time)<=adm.cens,avsurv$n.risk[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],0)
    
  } else {
    return(c(0,0,0,0,stopped=-1))
  }
  
  if (sum(data.sub$out)<5) stopped=1
  c(avsurv.est, avsurv.stderr, avsurv.num, dec.num, stopped) 
}#kmdec

# B. GND function
GND.calib = function(pred, tvar, out, cens.t, groups, adm.cens){
  output <- list()
  
  tvar.t=ifelse(tvar>adm.cens, adm.cens, tvar)
  out.t=ifelse(tvar>adm.cens, 0, out)
  
  datause=data.frame(pred=pred, tvar=tvar.t, out=out.t, count=1, cens.t=cens.t, dec=groups)
  numcat=length(unique(datause$dec))
  groups=sort(unique(datause$dec))
  
  kmtab=matrix(unlist(lapply(groups,kmdec,"dec",datain=datause, adm.cens)),ncol=5, byrow=TRUE)
  
  if (any(kmtab[,5] == -1)) stop("Stopped because at least one of the groups contains <2 events. Consider collapsing some groups.")
  else if (any(kmtab[,5] == 1)) warning("At least one of the groups contains < 5 events. GND can become unstable.\ 
(see Demler, Paynter, Cook 'Tests of Calibration and Goodness of Fit in the Survival Setting' DOI: 10.1002/sim.6428) \
Consider collapsing some groups to avoid this problem.")
  
  hltab=data.frame(group=kmtab[,4],
                   totaln=tapply(datause$count,datause$dec,sum),
                   censn=tapply(datause$cens.t,datause$dec,sum),
                   numevents=tapply(datause$out,datause$dec,sum),
                   expected=tapply(datause$pred,datause$dec,sum),
                   kmperc=1-kmtab[,1], 
                   kmvar=kmtab[,2]^2, 
                   kmnrisk=kmtab[,3],
                   expectedperc=tapply(datause$pred,datause$dec,mean))
  
  hltab$kmnum=hltab$kmperc*hltab$totaln
  hltab$GND_component=ifelse(hltab$kmvar==0, 0,(hltab$kmperc-hltab$expectedperc)^2/(hltab$kmvar))
  
  print(hltab[c(1,2,3,4,10,5,6,9,7,11)], digits=4)
  output[[1]] <- hltab[c(1,2,3,4,10,5,6,9,7,11)]
  output[[2]] <- c(df=numcat-1, chi2gw=sum(hltab$GND_component),pvalgw=1-pchisq(sum(hltab$GND_component),numcat-1))
  return(output)
}

# C. Quantile sorter
classifier <- function(risk,ncuts){
  cuts <- quantile(risk,probs=seq(0,1,1/ncuts))
  index <- rep(NA,length(risk))
  for (i in 1:(length(cuts)-1)){
    for (j in 1:length(risk)){
      index[j] <- ifelse(risk[j] >= cuts[i],i,index[j])}}
  return(index)
}

# D. Function to generate survival estimates per AF risk quantile
survivor <- function(data,risk_data,event,time,breakpoint){
  est <- rep(NA,times=length(unique(data[,get(risk_data)])))
  lower <- rep(NA,times=length(unique(data[,get(risk_data)])))
  upper <- rep(NA,times=length(unique(data[,get(risk_data)])))
  level_name <- rep(NA,times=length(unique(data[,get(risk_data)])))
  for (i in 1:length(unique(data[,get(risk_data)]))){
    subset <- data[get(risk_data)==unique(data[,get(risk_data)])[order(unique(data[,get(risk_data)]))][i]]
    if (nrow(subset[subset[,get(time)] >= breakpoint]) > 0){
      km <- survfit(Surv(subset[,get(time)],subset[,get(event)]) ~ 1, data=subset)
      time_index <- km$time - breakpoint
      end_time <- which(time_index == max(time_index[time_index <= 0]))
      est[i] <- 1-stepfun(km$time[1:end_time], c(1, km$surv[1:end_time]))(breakpoint)
      upper[i] <- 1-stepfun(km$time[1:end_time], c(1, km$lower[1:end_time]))(breakpoint)
      lower[i] <- 1-stepfun(km$time[1:end_time], c(1, km$upper[1:end_time]))(breakpoint)
      level_name[i] <- as.character(unique(data[,get(risk_data)])[order(unique(data[,get(risk_data)]))][i])
    }
    else {est[i] <- upper[i] <- lower[i] <- NA
    level_name[i] <- as.character(unique(data[,get(risk_data)])[order(unique(data[,get(risk_data)]))][i])}
  }
  return(data.frame(level=level_name,est=est,upper=upper,lower=lower))
}
