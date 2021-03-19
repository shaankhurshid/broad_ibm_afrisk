# Dependencies
library(data.table); library(survival); library(stringr)
source('~/broad_ibm_afrisk/subgroup_suite.R')

# Load data
load(file='/data/arrhythmia/skhurshid/p3po/charge_45.RData')

# Prev stroke
c3po_wide <- fread(file='/data/cvrepo/skhurshid/custom_wides/wide_with_mi_stroke_030121.csv')

numerics <- c('start_fu','ischemic_stroke_age')
for (j in numerics){set(c3po_wide,j=j,value=as.numeric(str_extract(c3po_wide[[j]],'\\d+')))}
c3po_wide[,prev_stroke := ifelse(c(!is.na(ischemic_stroke_age) & !is.na(start_fu) & (ischemic_stroke_age <= start_fu)),1,0)]
setkey(c3po_wide,id); setkey(charge_45,id)
charge_45[c3po_wide,prev_stroke := i.prev_stroke]

# FU in years
charge_45[,start_fu_yrs := start_fu/365.25]

## Define explore function for continuous variables (input must be data.table)
explore_age <- function(time,status,age_variable,min_age,max_age,
                        age_step,data,risk_score,path=getwd(),
                        pred_risk,threshold,calib_quantiles=10,censor.t,
                        make_plot=TRUE,all_pop=TRUE){
  i <- 1
  out <- list()
  for (age in seq(min_age,max_age-age_step,age_step)){
    
    if(age == max_age-age_step){
      subset <- data[c((get(age_variable) >= age) & (get(age_variable) <= age+age_step))] 
    } else { 
      subset <- data[c((get(age_variable) >= age) & (get(age_variable) < age+age_step))] }
    
    n_event <- nrow(subset[get(status)==1]); n_total <- nrow(subset)
    total_pt <- sum(subset[[time]]); event_ir <- n_event/total_pt
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
    
    # GND component
    gnd_obj <- GND.calib(pred=subset[[pred_risk]],tvar=subset[[time]],out=subset[[status]],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data="calib_groups",event=status,time=time,breakpoint=censor.t)
      obv <- incidence$est*100
      pred <- subset[,mean(get(pred_risk)),by='calib_groups']; setorder(pred,calib_groups)
      y_lim <- x_lim <- (max(obv,pred$V1,na.rm=T) %/% 5)*5+5
      pdf(file=paste0(path,'calib_',age,'.pdf'),height=3,width=3,pointsize=3)
      plot(pred$V1,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
      segments(-1,-1,101,101,lty=5)
      dev.off()}
    
    # Test char component
    if (tp == 0 | tn == 0 | test_pos == 0 | dz_pos == 0){
      ppv <- c(tp/test_pos,NA,NA); npv <- c(tn/test_neg,NA,NA); sens <- c(tp/dz_pos,NA,NA); spec <- c(tn/dz_neg,NA,NA)
    } else {
      ppv <- c(tp/test_pos,binom.test(tp,test_pos)$conf.int[1],binom.test(tp,test_pos)$conf.int[2])
      npv <- c(tn/test_neg,binom.test(tn,test_neg)$conf.int[1],binom.test(tn,test_neg)$conf.int[2])
      sens <- c(tp/dz_pos,binom.test(tp,dz_pos)$conf.int[1],binom.test(tp,dz_pos)$conf.int[2])
      spec <- c(tn/dz_neg,binom.test(tn,dz_neg)$conf.int[1],binom.test(tn,dz_neg)$conf.int[2])
    }
    
    # Writing output
    out[[i]] <- data.frame(matrix(ncol=30,nrow=0))
    out[[i]] <- as.numeric(c(paste0(age),n_event,n_total,total_pt,event_ir,
                             ci,as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                             as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3]),
                             relative_err,cum_err,gnd[2],gnd[3],
                             ppv,npv,sens,spec))
    names(out[[i]]) <- c('age','n_event','n_total','total_pt','event_ir',
                         'ci','ci_lower','ci_upper',
                         'c_stat','c_stat_lb','c_stat_ub',
                         'cal','cal_lb','cal_ub',
                         'rel_error_mean','abs_error_sum','gnd_chisq','gnd_p',
                         'ppv','ppv_lower','ppv_upper',
                         'npv','npv_lower','npv_upper',
                         'sens','sens_lower','sens_upper',
                         'spec','spec_lower','spec_upper')
    
    # Indicate progress
    print(paste0('Just finished model ',i,' out of ',((max_age-min_age)/age_step),'!'))
    i <- i+1
  }
  if (all_pop==TRUE){
    subset <- data
    n_event <- nrow(subset[get(status)==1]); n_total <- nrow(subset)
    total_pt <- sum(subset[[time]]); event_ir <- n_event/total_pt
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
    
    # GND component
    gnd_obj <- GND.calib(pred=subset[[pred_risk]],tvar=subset[[time]],out=subset[[status]],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data="calib_groups",event=status,time=time,breakpoint=censor.t)
      obv <- incidence$est*100
      pred <- subset[,mean(get(pred_risk)),by='calib_groups']; setorder(pred,calib_groups)
      y_lim <- x_lim <- (max(obv,pred$V1,na.rm=T) %/% 5)*5+5
      pdf(file=paste0(path,'calib_all.pdf'),height=3,width=3,pointsize=3)
      plot(pred$V1,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
      segments(-1,-1,101,101,lty=5)
      dev.off()}
    
    # Test char component
    if (tp == 0 | tn == 0 | test_pos == 0 | dz_pos == 0){
      ppv <- c(tp/test_pos,NA,NA); npv <- c(tn/test_neg,NA,NA); sens <- c(tp/dz_pos,NA,NA); spec <- c(tn/dz_neg,NA,NA)
    } else {
      ppv <- c(tp/test_pos,binom.test(tp,test_pos)$conf.int[1],binom.test(tp,test_pos)$conf.int[2])
      npv <- c(tn/test_neg,binom.test(tn,test_neg)$conf.int[1],binom.test(tn,test_neg)$conf.int[2])
      sens <- c(tp/dz_pos,binom.test(tp,dz_pos)$conf.int[1],binom.test(tp,dz_pos)$conf.int[2])
      spec <- c(tn/dz_neg,binom.test(tn,dz_neg)$conf.int[1],binom.test(tn,dz_neg)$conf.int[2])
    }
    
    # Writing output
    index <- length(out) + 1
    out[[index]] <- data.frame(matrix(ncol=30,nrow=0))
    out[[index]] <- c('All',as.numeric(c(n_event,n_total,total_pt,event_ir,
                                         ci,as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                                         as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3]),
                                         relative_err,cum_err,gnd[2],gnd[3],
                                         ppv,npv,sens,spec)))
    names(out[[index]]) <- c('age','n_event','n_total','total_pt','event_ir',
                             'ci','ci_lower','ci_upper',
                             'c_stat','c_stat_lb','c_stat_ub',
                             'cal','cal_lb','cal_ub',
                             'rel_error_mean','abs_error_sum','gnd_chisq','gnd_p',
                             'ppv','ppv_lower','ppv_upper',
                             'npv','npv_lower','npv_upper',
                             'sens','sens_lower','sens_upper',
                             'spec','spec_lower','spec_upper')
    print(paste0('Just finished all model!'))
  }
  return(data.frame(do.call(rbind,out)))
}

## Define explore function for categorical variables (input must be data.table)
explore_categorical <- function(time,status,variable,data,risk_score,path=getwd(),
                                pred_risk,threshold,calib_quantiles=10,censor.t,make_plot=TRUE,
                                all_pop=TRUE){
  i <- 1
  out <- list()
  for (var in unique(data[[variable]])){
    subset <- data[get(variable)==var]
    n_event <- nrow(subset[get(status)==1]); n_total <- nrow(subset)
    total_pt <- sum(subset[[time]]); event_ir <- n_event/total_pt
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
    
    # GND component
    gnd_obj <- GND.calib(pred=subset[[pred_risk]],tvar=subset[[time]],out=subset[[status]],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data="calib_groups",event=status,time=time,breakpoint=censor.t)
      obv <- incidence$est*100
      pred <- subset[,mean(get(pred_risk)),by='calib_groups']; setorder(pred,calib_groups)
      y_lim <- x_lim <- (max(obv,pred$V1) %/% 5)*5+5
      pdf(file=paste0(path,'calib_',variable,'_',var,'.pdf'),height=3,width=3,pointsize=3)
      plot(pred$V1,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
      segments(-1,-1,101,101,lty=5)
      dev.off()}
    
    if (tp == 0 | tn == 0 | test_pos == 0 | dz_pos == 0){
      ppv <- c(tp/test_pos,NA,NA); npv <- c(tn/test_neg,NA,NA); sens <- c(tp/dz_pos,NA,NA); spec <- c(tn/dz_neg,NA,NA)
    } else {
      ppv <- c(tp/test_pos,binom.test(tp,test_pos)$conf.int[1],binom.test(tp,test_pos)$conf.int[2])
      npv <- c(tn/test_neg,binom.test(tn,test_neg)$conf.int[1],binom.test(tn,test_neg)$conf.int[2])
      sens <- c(tp/dz_pos,binom.test(tp,dz_pos)$conf.int[1],binom.test(tp,dz_pos)$conf.int[2])
      spec <- c(tn/dz_neg,binom.test(tn,dz_neg)$conf.int[1],binom.test(tn,dz_neg)$conf.int[2])
    }
    out[[i]] <- data.frame(matrix(ncol=30,nrow=0))
    out[[i]] <- as.numeric(c(paste0(var),n_event,n_total,total_pt,event_ir,
                             ci,as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                             as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3]),
                             relative_err,cum_err,gnd[2],gnd[3],
                             ppv,npv,sens,spec))
    names(out[[i]]) <- c('stratum','n_event','n_total','total_pt','event_ir',
                         'ci','ci_lower','ci_upper',
                         'c_stat','c_stat_lb','c_stat_ub',
                         'cal','cal_lb','cal_ub',
                         'rel_error_mean','abs_error_sum','gnd_chisq','gnd_p',
                         'ppv','ppv_lower','ppv_upper',
                         'npv','npv_lower','npv_upper',
                         'sens','sens_lower','sens_upper',
                         'spec','spec_lower','spec_upper')
    print(paste0('Just finished model ',i,' out of ',length(unique(data[,get(variable)])),'!'))
    i <- i+1
  }
  if (all_pop==TRUE){
    subset <- data
    n_event <- nrow(subset[get(status)==1]); n_total <- nrow(subset)
    total_pt <- sum(subset[[time]]); event_ir <- n_event/total_pt
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
    
    # GND component
    gnd_obj <- GND.calib(pred=subset[[pred_risk]],tvar=subset[[time]],out=subset[[status]],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data="calib_groups",event=status,time=time,breakpoint=censor.t)
      obv <- incidence$est*100
      pred <- subset[,mean(get(pred_risk)),by='calib_groups']; setorder(pred,calib_groups)
      y_lim <- x_lim <- (max(obv,pred$V1) %/% 5)*5+5
      pdf(file=paste0(path,'calib_all.pdf'),height=3,width=3,pointsize=3)
      plot(pred$V1,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
      segments(-1,-1,101,101,lty=5)
      dev.off()}
    
    # Test char component
    if (tp == 0 | tn == 0 | test_pos == 0 | dz_pos == 0){
      ppv <- c(tp/test_pos,NA,NA); npv <- c(tn/test_neg,NA,NA); sens <- c(tp/dz_pos,NA,NA); spec <- c(tn/dz_neg,NA,NA)
    } else {
      ppv <- c(tp/test_pos,binom.test(tp,test_pos)$conf.int[1],binom.test(tp,test_pos)$conf.int[2])
      npv <- c(tn/test_neg,binom.test(tn,test_neg)$conf.int[1],binom.test(tn,test_neg)$conf.int[2])
      sens <- c(tp/dz_pos,binom.test(tp,dz_pos)$conf.int[1],binom.test(tp,dz_pos)$conf.int[2])
      spec <- c(tn/dz_neg,binom.test(tn,dz_neg)$conf.int[1],binom.test(tn,dz_neg)$conf.int[2])
    }
    
    # Writing output
    index <- length(out) + 1
    out[[index]] <- data.frame(matrix(ncol=30,nrow=0))
    out[[index]] <- c('All',as.numeric(c(n_event,n_total,total_pt,event_ir,
                                         ci,as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                                         as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3]),
                                         relative_err,cum_err,gnd[2],gnd[3],
                                         ppv,npv,sens,spec)))
    names(out[[index]]) <- c('age','n_event','n_total','total_pt','event_ir',
                             'ci','ci_lower','ci_upper',
                             'c_stat','c_stat_lb','c_stat_ub',
                             'cal','cal_lb','cal_ub',
                             'rel_error_mean','abs_error_sum','gnd_chisq','gnd_p',
                             'ppv','ppv_lower','ppv_upper',
                             'npv','npv_lower','npv_upper',
                             'sens','sens_lower','sens_upper',
                             'spec','spec_lower','spec_upper')
    print(paste0('Just finished all model!'))
  }
  return(data.frame(do.call(rbind,out)))
}

# Run explore functions
### Age
charge_age <- explore_age(time='af_5y.t',status='incd_af_5y',min_age=45,max_age=90,age_step=5,censor.t=5,
                          age_variable='start_fu_yrs',data=charge_45,path='/data/arrhythmia/skhurshid/heterogeneity/charge/',
                          risk_score='charge_startfu',pred_risk='charge_pred5',threshold=5,make_plot=TRUE,all_pop=FALSE)
write.csv(charge_age,file='/data/arrhythmia/skhurshid/heterogeneity/charge/charge_age.csv')

### Race
charge_race <- explore_categorical(time='af_5y.t',status='incd_af_5y',variable='race_white',
                                 data=charge_45,risk_score='charge_startfu',pred_risk='charge_pred5',censor.t=5,
                                 path='/data/arrhythmia/skhurshid/heterogeneity/charge/',threshold=5,make_plot=TRUE,all_pop=TRUE)
write.csv(charge_race,file='/data/arrhythmia/skhurshid/heterogeneity/charge/charge_race.csv')

### Sex
charge_sex <- explore_categorical(time='af_5y.t',status='incd_af_5y',variable='Dem.Gender.no_filter',
                                 data=charge_45,risk_score='charge_startfu',pred_risk='charge_pred5',censor.t=5,
                                 path='/data/arrhythmia/skhurshid/heterogeneity/charge/',threshold=5,make_plot=TRUE,all_pop=TRUE)
write.csv(charge_sex,file='/data/arrhythmia/skhurshid/heterogeneity/charge/charge_sex.csv')

### HF
charge_hf <- explore_categorical(time='af_5y.t',status='incd_af_5y',variable='prev_hf',
                                  data=charge_45,risk_score='charge_startfu',pred_risk='charge_pred5',censor.t=5,
                                  path='/data/arrhythmia/skhurshid/heterogeneity/charge/',threshold=5,make_plot=TRUE,all_pop=TRUE)
write.csv(charge_hf,file='/data/arrhythmia/skhurshid/heterogeneity/charge/charge_hf.csv')

### Stroke
charge_stroke <- explore_categorical(time='af_5y.t',status='incd_af_5y',variable='prev_stroke',
                                 data=charge_45,risk_score='charge_startfu',pred_risk='charge_pred5',censor.t=5,
                                 path='/data/arrhythmia/skhurshid/heterogeneity/charge/',threshold=5,make_plot=TRUE,all_pop=TRUE)
write.csv(charge_stroke,file='/data/arrhythmia/skhurshid/heterogeneity/charge/charge_stroke.csv')


