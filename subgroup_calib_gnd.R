# Dependencies
library(data.table); library(survival)

# Load data
load(file='/data/arrhythmia/skhurshid/ehr_af/vs_032120.RData')

# Script to explore coefficients in subgroups

# GND function
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

GND.calib = function(pred, tvar, out, cens.t, groups, adm.cens){
  
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
  
  c(df=numcat-1, chi2gw=sum(hltab$GND_component),pvalgw=1-pchisq(sum(hltab$GND_component),numcat-1))
}

# Quantile sorter
classifier <- function(risk,ncuts){
  cuts <- quantile(risk,probs=seq(0,1,1/ncuts))
  index <- rep(NA,length(risk))
  for (i in 1:(length(cuts)-1)){
    for (j in 1:length(risk)){
      index[j] <- ifelse(risk[j] >= cuts[i],i,index[j])}}
  return(index)
}

# Define explore function
explore_age <- function(time,status,age_variable,min_age,max_age,
                        age_step,data,risk_score,path=getwd(),
                        pred_risk,threshold,calib_quantiles=10,censor.t=5){
  i <- 1
  out <- list()
  for (age in seq(min_age,max_age,age_step)){
    if(age == max_age-age_step){
      subset <- data[c((get(age_variable) >= age) & (get(age_variable) < age+age_step))]
      n_af <- nrow(subset[get(status)==1]); n_total <- nrow(subset)
      mod <- coxph(Surv(subset[[time]],subset[[status]]) ~ subset[[risk_score]])
      km <- survfit(Surv(subset[[time]],subset[[status]]) ~ 1)
      ci <- c((1-km$surv[length(km$surv)])*100,
              (1-km$upper[length(km$upper)])*100,
              (1-km$lower[length(km$lower)])*100)
      tp <- nrow(subset[c((get(pred_risk) >= threshold) & get(status)==1)])
      tn <- nrow(subset[c((get(pred_risk) < threshold) & get(status)==0)])
      test_pos <- nrow(subset[get(pred_risk) >= threshold]); test_neg <- nrow(subset[get(pred_risk) < threshold])
      dz_pos <- nrow(subset[get(status)==1]); dz_neg <- nrow(subset[get(status)==0])
      subset[,':='(calib_groups = classifier(risk=get(pred_risk),ncuts=calib_quantiles),
                   censored = ifelse(get(status)==0,1,0))]
      gnd <- GND.calib(pred=subset[[pred_risk]],tvar=subset[[time]],out=subset[[status]],groups=subset$calib_groups,
                       cens.t=subset$censored,adm.cens = censor.t)
      
      obv <- subset[,lapply(.SD,mean),by='calib_groups',.SDcols=status][order(calib_groups),get(status)]*100
      pred <- subset[,mean(get(pred_risk)),by='calib_groups'][order(calib_groups)]
      
      pdf(file=paste0(path,'calib_',age,'.pdf'),height=3,width=3,pointsize=3)
      plot(pred$V1*100,obv,xlab='Predicted',ylab='Observed',xlim=c(0,30),ylim=c(0,30),pch=19)
      segments(-1,-1,101,101,lty=5)
      dev.off()
      
      if (tp == 0 | tn == 0 | test_pos == 0 | dz_pos == 0){
        ppv <- c(tp/test_pos,NA,NA); npv <- c(tn/test_neg,NA,NA); sens <- c(tp/dz_pos,NA,NA); spec <- c(tn/dz_neg,NA,NA)
      } else {
        ppv <- c(tp/test_pos,binom.test(tp,test_pos)$conf.int[1],binom.test(tp,test_pos)$conf.int[2])
        npv <- c(tn/test_neg,binom.test(tn,test_neg)$conf.int[1],binom.test(tn,test_neg)$conf.int[2])
        sens <- c(tp/dz_pos,binom.test(tp,dz_pos)$conf.int[1],binom.test(tp,dz_pos)$conf.int[2])
        spec <- c(tn/dz_neg,binom.test(tn,dz_neg)$conf.int[1],binom.test(tn,dz_neg)$conf.int[2])
      }
      out[[i]] <- data.frame(matrix(ncol=26,nrow=0))
      out[[i]] <- as.numeric(c(paste0(age),n_af,n_total,ci,
                               as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                               as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3]),
                               gnd[2],gnd[3],
                               ppv,npv,sens,spec))
      names(out[[i]]) <- c('age','n_af','n_total',
                           'ci','ci_lower','ci_upper',
                           'c_stat','c_stat_lb','c_stat_ub',
                           'cal','cal_lb','cal_ub',
                           'gnd_chisq','gnd_p',
                           'ppv','ppv_lower','ppv_upper',
                           'npv','npv_lower','npv_upper',
                           'sens','sens_lower','sens_upper',
                           'spec','spec_lower','spec_upper')
      print(paste0('Just finished model ',i,' out of ',((max_age-min_age)/age_step),'!'))
      break
    } else {
      subset <- data[c((get(age_variable) >= age) & (get(age_variable) < age+age_step))]
      n_af <- nrow(subset[get(status)==1]); n_total <- nrow(subset)
      mod <- coxph(Surv(subset[[time]],subset[[status]]) ~ subset[[risk_score]])
      km <- survfit(Surv(subset[[time]],subset[[status]]) ~ 1)
      ci <- c((1-km$surv[length(km$surv)])*100,
              (1-km$upper[length(km$upper)])*100,
              (1-km$lower[length(km$lower)])*100)
      tp <- nrow(subset[c((get(pred_risk) >= threshold) & get(status)==1)])
      tn <- nrow(subset[c((get(pred_risk) < threshold) & get(status)==0)])
      test_pos <- nrow(subset[get(pred_risk) >= threshold]); test_neg <- nrow(subset[get(pred_risk) < threshold])
      dz_pos <- nrow(subset[get(status)==1]); dz_neg <- nrow(subset[get(status)==0])
      subset[,':='(calib_groups = classifier(risk=get(pred_risk),ncuts=calib_quantiles),
                   censored = ifelse(get(status)==0,1,0))]
      gnd <- GND.calib(pred=subset[[pred_risk]],tvar=subset[[time]],out=subset[[status]],groups=subset$calib_groups,
                       cens.t=subset$censored,adm.cens = censor.t)
      
      obv <- subset[,lapply(.SD,mean),by='calib_groups',.SDcols=status][order(calib_groups),get(status)]*100
      pred <- subset[,mean(get(pred_risk)),by='calib_groups'][order(calib_groups)]
      
      pdf(file=paste0(path,'calib_',age,'.pdf'),height=3,width=3,pointsize=3)
      plot(pred$V1*100,obv,xlab='Predicted',ylab='Observed',xlim=c(0,30),ylim=c(0,30),pch=19)
      segments(-1,-1,101,101,lty=5)
      dev.off()
      
      if (tp == 0 | tn == 0 | test_pos == 0 | dz_pos == 0){
        ppv <- c(tp/test_pos,NA,NA); npv <- c(tn/test_neg,NA,NA); sens <- c(tp/dz_pos,NA,NA); spec <- c(tn/dz_neg,NA,NA)
      } else {
        ppv <- c(tp/test_pos,binom.test(tp,test_pos)$conf.int[1],binom.test(tp,test_pos)$conf.int[2])
        npv <- c(tn/test_neg,binom.test(tn,test_neg)$conf.int[1],binom.test(tn,test_neg)$conf.int[2])
        sens <- c(tp/dz_pos,binom.test(tp,dz_pos)$conf.int[1],binom.test(tp,dz_pos)$conf.int[2])
        spec <- c(tn/dz_neg,binom.test(tn,dz_neg)$conf.int[1],binom.test(tn,dz_neg)$conf.int[2])
      }
      out[[i]] <- data.frame(matrix(ncol=26,nrow=0))
      out[[i]] <- as.numeric(c(paste0(age),n_af,n_total,ci,
                               as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                               as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3]),
                               gnd[2],gnd[3],
                               ppv,npv,sens,spec))
      names(out[[i]]) <- c('age','n_af','n_total',
                           'ci','ci_lower','ci_upper',
                           'c_stat','c_stat_lb','c_stat_ub',
                           'cal','cal_lb','cal_ub',
                           'gnd_chisq','gnd_p',
                           'ppv','ppv_lower','ppv_upper',
                           'npv','npv_lower','npv_upper',
                           'sens','sens_lower','sens_upper',
                           'spec','spec_lower','spec_upper')
      print(paste0('Just finished model ',i,' out of ',((max_age-min_age)/age_step),'!'))
      i <- i+1
    }
  }
  return(data.frame(do.call(rbind,out)))
}

# Define explore function
explore_categorical <- function(time,status,variable,data,risk_score,path=getwd(),
                                pred_risk,threshold,calib_quantiles=10,censor.t=5){
  i <- 1
  out <- list()
  for (var in unique(data[[variable]])){
    subset <- data[get(variable)==var]
    n_af <- nrow(subset[get(status)==1]); n_total <- nrow(subset)
    mod <- coxph(Surv(subset[[time]],subset[[status]]) ~ subset[[risk_score]])
    km <- survfit(Surv(subset[[time]],subset[[status]]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    tp <- nrow(subset[c((get(pred_risk) >= threshold) & get(status)==1)])
    tn <- nrow(subset[c((get(pred_risk) < threshold) & get(status)==0)])
    test_pos <- nrow(subset[get(pred_risk) >= threshold]); test_neg <- nrow(subset[get(pred_risk) < threshold])
    dz_pos <- nrow(subset[get(status)==1]); dz_neg <- nrow(subset[get(status)==0])
    subset[,':='(calib_groups = classifier(risk=get(pred_risk),ncuts=calib_quantiles),
                 censored = ifelse(get(status)==0,1,0))]
    gnd <- GND.calib(pred=subset[[pred_risk]],tvar=subset[[time]],out=subset[[status]],groups=subset$calib_groups,
                     cens.t=subset$censored,adm.cens = censor.t)
    
    obv <- subset[,lapply(.SD,mean),by='calib_groups',.SDcols=status][order(calib_groups),get(status)]*100
    pred <- subset[,mean(get(pred_risk)),by='calib_groups'][order(calib_groups)]
    
    pdf(file=paste0(path,'calib_',var,'.pdf'),height=3,width=3,pointsize=3)
    plot(pred$V1*100,obv,xlab='Predicted',ylab='Observed',xlim=c(0,30),ylim=c(0,30),pch=19)
    segments(-1,-1,101,101,lty=5)
    dev.off()
    
    if (tp == 0 | tn == 0 | test_pos == 0 | dz_pos == 0){
      ppv <- c(tp/test_pos,NA,NA); npv <- c(tn/test_neg,NA,NA); sens <- c(tp/dz_pos,NA,NA); spec <- c(tn/dz_neg,NA,NA)
    } else {
      ppv <- c(tp/test_pos,binom.test(tp,test_pos)$conf.int[1],binom.test(tp,test_pos)$conf.int[2])
      npv <- c(tn/test_neg,binom.test(tn,test_neg)$conf.int[1],binom.test(tn,test_neg)$conf.int[2])
      sens <- c(tp/dz_pos,binom.test(tp,dz_pos)$conf.int[1],binom.test(tp,dz_pos)$conf.int[2])
      spec <- c(tn/dz_neg,binom.test(tn,dz_neg)$conf.int[1],binom.test(tn,dz_neg)$conf.int[2])
    }
    out[[i]] <- data.frame(matrix(ncol=26,nrow=0))
    out[[i]] <- as.numeric(c(paste0(var),n_af,n_total,ci,
                             as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                             as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3]),
                             gnd[2],gnd[3],
                             ppv,npv,sens,spec))
    names(out[[i]]) <- c('category','n_af','n_total',
                         'ci','ci_lower','ci_upper',
                         'c_stat','c_stat_lb','c_stat_ub',
                         'cal','cal_lb','cal_ub',
                         'ppv','ppv_lower','ppv_upper',
                         'npv','npv_lower','npv_upper',
                         'sens','sens_lower','sens_upper',
                         'spec','spec_lower','spec_upper')
    print(paste0('Just finished model ',i,' out of ',length(unique(data[,get(variable)])),'!'))
    i <- i+1
  }
  return(data.frame(do.call(rbind,out)))
}

# Run explore functions
output_age <- explore_age(time='af_5y_sal.t',status='af_5y_sal',min_age=45,max_age=95,age_step=5,
                          age_variable='start_fu_age',data=vs,path='/data/arrhythmia/skhurshid/broad_ibm_afrisk/',
                          risk_score='score',pred_risk='pred5',threshold=0.05)

output_sex <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='Gender',
                                  data=vs,risk_score='score',pred_risk='pred5',
                                  path='/data/arrhythmia/skhurshid/broad_ibm_afrisk/',threshold=0.05)
