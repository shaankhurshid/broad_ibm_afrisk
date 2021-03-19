# Dependencies
library(survival); library(data.table)

# Load data
#load(file='/data/arrhythmia/skhurshid/ehr_af/vs_032120.RData')
vs <- fread('/data/arrhythmia/skhurshid/heterogeneity/vs_439k.csv')
setDT(vs)

# Source helper functions
source('subgroup_suite_df.R')

names(vs)[names(vs) == "is_AFib"] <- "af_5y_sal"
names(vs)[names(vs) == "num_days_to_AFib_or_to_cencor"] <- "af_5y_sal.t"
names(vs)[names(vs) == "CHARGE_AF"] <- "chargeaf"
names(vs)[names(vs) == "age_at_index_date"] <- "start_fu_age"
names(vs)[names(vs) == "Gender"] <- "is_male"

# Script to explore coefficients in subgroups

## Define explore function for continuous variables (input must be data.table)
explore_age <- function(time,status,age_variable,min_age,max_age,
                        age_step,data,risk_score,path=getwd(),
                        pred_risk,threshold,calib_quantiles=10,censor.t=5,
                        make_plot=TRUE,all_pop=TRUE){
  i <- 1
  out <- list()
  for (age in seq(min_age,max_age-age_step,age_step)){
    
    if(age == max_age-age_step){
      subset <- data[c((data[,age_variable] >= age) & (data[,age_variable] <= age+age_step)),] 
    } else { 
      subset <- data[c((data[,age_variable] >= age) & (data[,age_variable] < age+age_step)),]}
    
    n_event <- nrow(subset[subset[,status]==1,]); n_total <- nrow(subset)
    total_pt <- sum(subset[,time]); event_ir <- n_event/total_pt
    mod <- coxph(Surv(subset[,time],subset[,status]) ~ subset[,risk_score])
    km <- survfit(Surv(subset[,time],subset[,status]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    tp <- nrow(subset[c((subset[,pred_risk] >= threshold) & subset[,status]==1),])
    tn <- nrow(subset[c((subset[,pred_risk] < threshold) & subset[,status]==0),])
    test_pos <- nrow(subset[subset[,pred_risk] >= threshold,]); test_neg <- nrow(subset[subset[,pred_risk] < threshold,])
    dz_pos <- nrow(subset[subset[,status]==1,]); dz_neg <- nrow(subset[subset[,status]==0,])
    subset$calib_groups <- classifier(risk=subset[,pred_risk],ncuts=calib_quantiles)
    subset$censored <- ifelse(subset[,status]==0,1,0)
    
    # GND component
    gnd_obj <- GND.calib(pred=subset[,pred_risk],tvar=subset[,time],out=subset[,status],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data="calib_groups",event=status,time=time,breakpoint=5)
      obv <- incidence$est*100
      pred <- unlist(lapply(split(subset[,pred_risk],subset$calib_groups),mean,na.rm=TRUE))
      y_lim <- x_lim <- (max(obv,pred*100,na.rm=T) %/% 5)*5+5
      pdf(file=paste0(path,'calib_',age,'.pdf'),height=3,width=3,pointsize=3)
      plot(pred*100,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
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
    n_event <- nrow(subset[subset[,status]==1,]); n_total <- nrow(subset)
    total_pt <- sum(subset[,time]); event_ir <- n_event/total_pt
    mod <- coxph(Surv(subset[,time],subset[,status]) ~ subset[,risk_score])
    km <- survfit(Surv(subset[,time],subset[,status]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    tp <- nrow(subset[c((subset[,pred_risk] >= threshold) & subset[,status]==1),])
    tn <- nrow(subset[c((subset[,pred_risk] < threshold) & subset[,status]==0),])
    test_pos <- nrow(subset[subset[,pred_risk] >= threshold,]); test_neg <- nrow(subset[subset[,pred_risk] < threshold,])
    dz_pos <- nrow(subset[subset[,status]==1,]); dz_neg <- nrow(subset[subset[,status]==0,])
    subset$calib_groups <- classifier(risk=subset[,pred_risk],ncuts=calib_quantiles)
    subset$censored <- ifelse(subset[,status]==0,1,0)
    
    # GND component
    gnd_obj <- GND.calib(pred=subset[,pred_risk],tvar=subset[,time],out=subset[,status],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data="calib_groups",event=status,time=time,breakpoint=5)
      obv <- incidence$est*100
      pred <- unlist(lapply(split(subset[,pred_risk],subset$calib_groups),mean,na.rm=TRUE))
      y_lim <- x_lim <- (max(obv,pred*100,na.rm=T) %/% 5)*5+5
      pdf(file=paste0(path,'calib_',age,'.pdf'),height=3,width=3,pointsize=3)
      plot(pred*100,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
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
                                pred_risk,threshold,calib_quantiles=10,censor.t=5,make_plot=TRUE,
                                all_pop=TRUE){
  i <- 1
  out <- list()
  for (var in unique(data[,variable])){
    subset <- data[data[,variable]==var,]
    n_event <- nrow(subset[subset[,status]==1,]); n_total <- nrow(subset)
    total_pt <- sum(subset[,time]); event_ir <- n_event/total_pt
    mod <- coxph(Surv(subset[,time],subset[,status]) ~ subset[,risk_score])
    km <- survfit(Surv(subset[,time],subset[,status]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    tp <- nrow(subset[c((subset[,pred_risk] >= threshold) & subset[,status]==1),])
    tn <- nrow(subset[c((subset[,pred_risk] < threshold) & subset[,status]==0),])
    test_pos <- nrow(subset[subset[,pred_risk] >= threshold,]); test_neg <- nrow(subset[subset[,pred_risk] < threshold,])
    dz_pos <- nrow(subset[subset[,status]==1,]); dz_neg <- nrow(subset[subset[,status]==0,])
    subset$calib_groups <- classifier(risk=subset[,pred_risk],ncuts=calib_quantiles)
    subset$censored <- ifelse(subset[,status]==0,1,0)
    
    # GND component
    gnd_obj <- GND.calib(pred=subset[,pred_risk],tvar=subset[,time],out=subset[,status],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data="calib_groups",event=status,time=time,breakpoint=5)
      obv <- incidence$est*100
      pred <- unlist(lapply(split(subset[,pred_risk],subset$calib_groups),mean,na.rm=TRUE))
      y_lim <- x_lim <- (max(obv,pred*100,na.rm=T) %/% 5)*5+5
      pdf(file=paste0(path,'calib_',var,'.pdf'),height=3,width=3,pointsize=3)
      plot(pred*100,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
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
    
    # Writing component
    out[[i]] <- data.frame(matrix(ncol=30,nrow=0))
    out[[i]] <- c(paste0(var),as.numeric(c(n_event,n_total,total_pt,event_ir,
                                           ci,as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                                           as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3]),
                                           relative_err,cum_err,gnd[2],gnd[3],
                                           ppv,npv,sens,spec)))
    names(out[[i]]) <- c('stratum','n_event','n_total','total_pt','event_ir',
                         'ci','ci_lower','ci_upper',
                         'c_stat','c_stat_lb','c_stat_ub',
                         'cal','cal_lb','cal_ub',
                         'rel_error_mean','abs_error_sum','gnd_chisq','gnd_p',
                         'ppv','ppv_lower','ppv_upper',
                         'npv','npv_lower','npv_upper',
                         'sens','sens_lower','sens_upper',
                         'spec','spec_lower','spec_upper')
    print(paste0('Just finished model ',i,' out of ',length(unique(data[,variable])),'!'))
    i <- i+1
  }
  
  if (all_pop==TRUE){
    subset <- data
    n_event <- nrow(subset[subset[,status]==1,]); n_total <- nrow(subset)
    total_pt <- sum(subset[,time]); event_ir <- n_event/total_pt
    mod <- coxph(Surv(subset[,time],subset[,status]) ~ subset[,risk_score])
    km <- survfit(Surv(subset[,time],subset[,status]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    tp <- nrow(subset[c((subset[,pred_risk] >= threshold) & subset[,status]==1),])
    tn <- nrow(subset[c((subset[,pred_risk] < threshold) & subset[,status]==0),])
    test_pos <- nrow(subset[subset[,pred_risk] >= threshold,]); test_neg <- nrow(subset[subset[,pred_risk] < threshold,])
    dz_pos <- nrow(subset[subset[,status]==1,]); dz_neg <- nrow(subset[subset[,status]==0,])
    subset$calib_groups <- classifier(risk=subset[,pred_risk],ncuts=calib_quantiles)
    subset$censored <- ifelse(subset[,status]==0,1,0)
    
    # GND component
    gnd_obj <- GND.calib(pred=subset[,pred_risk],tvar=subset[,time],out=subset[,status],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data="calib_groups",event=status,time=time,breakpoint=5)
      obv <- incidence$est*100
      pred <- unlist(lapply(split(subset[,pred_risk],subset$calib_groups),mean,na.rm=TRUE))
      y_lim <- x_lim <- (max(obv,pred*100,na.rm=T) %/% 5)*5+5
      pdf(file=paste0(path,'calib_',var,'.pdf'),height=3,width=3,pointsize=3)
      plot(pred*100,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
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


#Converting time to event from days to years.
vs$af_5y_sal.t <- vs$af_5y_sal.t / 365
vs$af_5y_sal.t <- as.numeric(vs$af_5y_sal.t)

vs[,score_avgbeta_chargeaf := mean(chargeaf)]
res_chargeaf <- coxph(Surv(af_5y_sal.t,af_5y_sal) ~ chargeaf, data = vs); summary(res_chargeaf)
km_chargeaf <- survfit(res_chargeaf, data=data.frame(x1 = mean(chargeaf)),type="kaplan-meier")
s0_chargeaf <- summary(km_chargeaf, times = c(5))$surv; s0_chargeaf
vs[,charge.pred5 := (1-(s0_chargeaf)^exp(chargeaf - (score_avgbeta_chargeaf)))]; summary(vs$charge.pred5)

#head(vs$charge.pred5)

setDF(vs)
# Run explore functions
output_age <- explore_age(time='af_5y_sal.t',status='af_5y_sal',min_age=45,max_age=55,age_step=5,
                          age_variable='start_fu_age',data=vs,path='/data/arrhythmia/skhurshid/heterogeneity/test/',
                          risk_score='chargeaf',pred_risk='charge.pred5',threshold=0.05,make_plot=TRUE,all_pop=TRUE)
write.csv(output_age,file='//charge_output_age.csv')

output_sex <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='Gender',
                                  data=vs,risk_score='chargeaf',pred_risk='charge.pred5',
                                  path='/',threshold=0.05,make_plot=TRUE,all_pop=TRUE)
write.csv(output_sex,file='/charge_output_sex.csv')

output_hf <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='heartFailure_fin_prevAtstartFu',
                                 data=vs,risk_score='chargeaf',pred_risk='charge.pred5',
                                 path='/',threshold=0.05,make_plot=TRUE,all_pop=TRUE)
write.csv(output_sex,file='/charge_output_hf.csv')

output_stroke <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='cvaTia_fin_prevAtstartFu',
                                     data=vs,risk_score='chargeaf',pred_risk='charge.pred5',
                                     path='/',threshold=0.05,make_plot=TRUE,all_pop=TRUE)
write.csv(output_sex,file='/charge_output_stroke.csv')

output_race <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='race_binary',
                                   data=vs,risk_score='chargeaf',pred_risk='charge.pred5',
                                   path='/',threshold=0.05,make_plot=TRUE,all_pop=TRUE)
write.csv(output_race,file='/charge_output_race.csv')
