# Dependencies
library(data.table); library(survival); library(stringr)
source('~/broad_ibm_afrisk/subgroup_suite.R')

# Load data
load(file='/data/arrhythmia/skhurshid/heterogeneity/pce_40.RData')

# Prev HF
numerics <- c('start_fu','inpatient_hf_age')
for (j in numerics){set(pce_40,j=j,value=as.numeric(str_extract(pce_40[[j]],'\\d+')))}
pce_40[,prev_hf := ifelse(c(!is.na(inpatient_hf_age) & !is.na(start_fu) & (inpatient_hf_age <= start_fu)),1,0)]

# Subset data for convenience
pce_wf <- pce_40[c(Dem.Gender.no_filter=='Female' & race_black==0)]
pce_bf <- pce_40[c(Dem.Gender.no_filter=='Female' & race_black==1)]
pce_wm <- pce_40[c(Dem.Gender.no_filter=='Male' & race_black==0)]
pce_bm <- pce_40[c(Dem.Gender.no_filter=='Male' & race_black==1)]

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
age_wf <- explore_age(time='mi_stroke_10y.t',status='incd_mi_stroke_10y',min_age=40,max_age=80,age_step=10,censor.t=10,
                          age_variable='start_fu_yrs',data=pce_wf,path='/data/arrhythmia/skhurshid/heterogeneity/pce_wf/',
                          risk_score='pce',pred_risk='pce_risk',threshold=7.5,make_plot=TRUE,all_pop=TRUE)
write.csv(age_wf,file='/data/arrhythmia/skhurshid/heterogeneity/pce_wf/pce_age_wf.csv')

age_bf <- explore_age(time='mi_stroke_10y.t',status='incd_mi_stroke_10y',min_age=40,max_age=80,age_step=10,censor.t=10,
                      age_variable='start_fu_yrs',data=pce_bf,path='/data/arrhythmia/skhurshid/heterogeneity/pce_bf/',
                      risk_score='pce',pred_risk='pce_risk',threshold=7.5,make_plot=TRUE,all_pop=TRUE)
write.csv(age_bf,file='/data/arrhythmia/skhurshid/heterogeneity/pce_bf/pce_age_bf.csv')

age_wm <- explore_age(time='mi_stroke_10y.t',status='incd_mi_stroke_10y',min_age=40,max_age=80,age_step=10,censor.t=10,
                      age_variable='start_fu_yrs',data=pce_wm,path='/data/arrhythmia/skhurshid/heterogeneity/pce_wm/',
                      risk_score='pce',pred_risk='pce_risk',threshold=7.5,make_plot=TRUE,all_pop=TRUE)
write.csv(age_wm,file='/data/arrhythmia/skhurshid/heterogeneity/pce_wm/pce_age_wm.csv')

age_bm <- explore_age(time='mi_stroke_10y.t',status='incd_mi_stroke_10y',min_age=40,max_age=80,age_step=10,censor.t=10,
                      age_variable='start_fu_yrs',data=pce_bm,path='/data/arrhythmia/skhurshid/heterogeneity/pce_bm/',
                      risk_score='pce',pred_risk='pce_risk',threshold=7.5,make_plot=TRUE,all_pop=TRUE,calib_quantiles=5)
write.csv(age_bm,file='/data/arrhythmia/skhurshid/heterogeneity/pce_bm/pce_age_bm.csv')

### HF
hf_wf <- explore_categorical(time='mi_stroke_10y.t',status='incd_mi_stroke_10y',variable='prev_hf',
                                  data=pce_wf,risk_score='pce',pred_risk='pce_risk',censor.t=10,
                                  path='/data/arrhythmia/skhurshid/heterogeneity/pce_wf/',threshold=7.5,make_plot=TRUE,all_pop=TRUE)
write.csv(hf_wf,file='/data/arrhythmia/skhurshid/heterogeneity/pce_wf/pce_hf_wf.csv')

hf_bf <- explore_categorical(time='mi_stroke_10y.t',status='incd_mi_stroke_10y',variable='prev_hf',
                             data=pce_bf,risk_score='pce',pred_risk='pce_risk',censor.t=10, calib_quantiles=5,
                             path='/data/arrhythmia/skhurshid/heterogeneity/pce_bf/',threshold=7.5,make_plot=TRUE,all_pop=TRUE)
write.csv(hf_bf,file='/data/arrhythmia/skhurshid/heterogeneity/pce_bf/pce_hf_bf.csv')

hf_wm <- explore_categorical(time='mi_stroke_10y.t',status='incd_mi_stroke_10y',variable='prev_hf',
                             data=pce_wm,risk_score='pce',pred_risk='pce_risk',censor.t=10,
                             path='/data/arrhythmia/skhurshid/heterogeneity/pce_wm/',threshold=7.5,make_plot=TRUE,all_pop=TRUE)
write.csv(hf_wm,file='/data/arrhythmia/skhurshid/heterogeneity/pce_wm/pce_hf_wm.csv')

hf_bm <- explore_categorical(time='mi_stroke_10y.t',status='incd_mi_stroke_10y',variable='prev_hf',
                             data=pce_bm,risk_score='pce',pred_risk='pce_risk',censor.t=10, calib_quantiles=5,
                             path='/data/arrhythmia/skhurshid/heterogeneity/pce_bm/',threshold=7.5,make_plot=TRUE,all_pop=TRUE)
write.csv(hf_bm,file='/data/arrhythmia/skhurshid/heterogeneity/pce_bm/pce_hf_bm.csv')
