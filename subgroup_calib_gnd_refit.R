# Dependencies
library(data.table); library(survival)

# Source helper functions
source('subgroup_suite.R')

# Load data
load(file='/data/arrhythmia/skhurshid/ehr_af/vs_032120.RData')
setDT(vs)

## Define explore function for continuous variables (input must be data.table)
explore_continuous_recal <- function(time,status,age_variable,min_age,max_age,variables_list,
                            age_step,dataset,path=getwd(),threshold,
                            calib_quantiles=10,censor.t=5,
                            make_plot=TRUE,all_pop=TRUE){
  i <- 1
  out <- list()
  for (age in seq(min_age,max_age-age_step,age_step)){
    
    if(age == max_age-age_step){
      subset <- dataset[c((get(age_variable) >= age) & (get(age_variable) <= age+age_step))] 
    } else { 
      subset <- dataset[c((get(age_variable) >= age) & (get(age_variable) < age+age_step))] }
      n_event <- nrow(subset[get(status)==1]); n_total <- nrow(subset)
      
      total_pt <- sum(subset[[time]]); event_ir <- n_event/total_pt
      km <- survfit(Surv(subset[[time]],subset[[status]]) ~ 1)
      ci <- c((1-km$surv[length(km$surv)])*100,
              (1-km$upper[length(km$upper)])*100,
              (1-km$lower[length(km$lower)])*100)
      
      # Refit model
      new_formula <- formula(paste0('Surv(subset[,get(time)],subset[,get(status)]) ~ ',paste0(variables_list,collapse=' + ')))
      mod <- coxph(new_formula,data=subset)
      subset[,lp := predict(mod,type='lp')]
      
      # Calculate predicted risk with recalibration
      avg_beta <- mean(subset$lp)
      res <- coxph(Surv(subset[,time],subset[,status]) ~ lp, data=subset)
      km <- survfit(res, data=data.frame(x1=mean(lp)),type="kaplan-meier")
      s0 <- summary(km, times=c(censor.t))$surv
      subset[,pred_risk := (1-(s0)^exp(lp - (avg_beta)))]

      # Sens/spec metrics
      tp <- nrow(subset[c(pred_risk >= threshold & get(status)==1)])
      tn <- nrow(subset[c(pred_risk < threshold & get(status)==0)])
      test_pos <- nrow(subset[pred_risk >= threshold]); test_neg <- nrow(subset[pred_risk < threshold])
      dz_pos <- nrow(subset[get(status)==1]); dz_neg <- nrow(subset[get(status)==0])
      
      subset[,':='(calib_groups = classifier(risk=subset$pred_risk,ncuts=calib_quantiles),
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
        incidence <- survivor(data=subset,risk_data="calib_groups",event=event,time=time,breakpoint=5)
        obv <- incidence$est*100
        pred <- subset[,mean(pred_risk),by='calib_groups']; setorder(pred,calib_groups)
        y_lim <- x_lim <- (max(obv,pred$V1*100,na.rm=T) %/% 5)*5+5
        pdf(file=paste0(path,'recalib_',age,'.pdf'),height=3,width=3,pointsize=3)
        plot(pred$V1*100,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
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
      out[[i]] <- data.frame(matrix(ncol=32,nrow=0))
      out[[i]] <- as.numeric(c(paste0(age),n_event,n_total,total_pt,event_ir,
                               ci,as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                               s0,avg_beta,
                               as.numeric(res$coefficients[1]),as.numeric(res$coefficients[1]-1.96*summary(res)$coefficients[3]),as.numeric(res$coefficients[1]+1.96*summary(res)$coefficients[3]),
                               relative_err,cum_err,gnd[2],gnd[3],
                               ppv,npv,sens,spec))
      names(out[[i]]) <- c('age','n_event','n_total','total_pt','event_ir',
                           'ci','ci_lower','ci_upper',
                           'c_stat','c_stat_lb','c_stat_ub',
                           's0','avg_beta',
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
    subset <- dataset
    n_event <- nrow(subset[get(status)==1]); n_total <- nrow(subset)
    total_pt <- sum(subset[[time]]); event_ir <- n_event/total_pt
    mod <- coxph(Surv(subset[[time]],subset[[status]]) ~ subset[[risk_score]])
    km <- survfit(Surv(subset[[time]],subset[[status]]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    
    # Refit model
    new_formula <- formula(paste0('Surv(subset[,get(time)],subset[,get(status)]) ~ ',paste0(variables_list,collapse=' + ')))
    mod <- coxph(new_formula,data=subset)
    subset[,lp := predict(mod,type='lp')]
    
    # Calculate predicted risk with recalibration
    avg_beta <- mean(subset$lp)
    res <- coxph(Surv(subset[,get(time)],subset[,get(status)]) ~ lp, data=subset)
    km <- survfit(res, data=data.frame(x1=mean(lp)),type="kaplan-meier")
    s0 <- summary(km, times=c(5))$surv
    subset[,pred_risk := (1-(s0)^exp(lp - (avg_beta)))]
    
    # Sens/spec tests
    tp <- nrow(subset[c(pred_risk >= threshold & get(status)==1)])
    tn <- nrow(subset[c(pred_risk < threshold & get(status)==0)])
    test_pos <- nrow(subset[pred_risk >= threshold]); test_neg <- nrow(subset[pred_risk < threshold])
    dz_pos <- nrow(subset[get(status)==1]); dz_neg <- nrow(subset[get(status)==0])
    subset[,':='(calib_groups = classifier(risk=subset$pred_risk,ncuts=calib_quantiles),
                 censored = ifelse(get(status)==0,1,0))]
    
    # GND component
    gnd_obj <- GND.calib(pred=subset[[pred_risk]],tvar=subset[[time]],out=subset[[status]],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      incidence <- survivor(data=subset,risk_data="calib_groups",event=event,time=time,breakpoint=5)
      obv <- incidence$est*100
      pred <- subset[,mean(pred_risk),by='calib_groups']; setorder(pred,calib_groups)
      y_lim <- x_lim <- (max(obv,pred$V1*100,na.rm=T) %/% 5)*5+5
      pdf(file=paste0(path,'recalib_all.pdf'),height=3,width=3,pointsize=3)
      plot(pred$V1*100,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
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
    out[[index]] <- data.frame(matrix(ncol=32,nrow=0))
    out[[index]] <- c('All',as.numeric(c(n_event,n_total,total_pt,event_ir,
                             ci,as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                             s0,avg_beta,
                             as.numeric(res$coefficients[1]),as.numeric(res$coefficients[1]-1.96*summary(res)$coefficients[3]),as.numeric(res$coefficients[1]+1.96*summary(res)$coefficients[3]),
                             relative_err,cum_err,gnd[2],gnd[3],
                             ppv,npv,sens,spec)))
    names(out[[index]]) <- c('age','n_event','n_total','total_pt','event_ir',
                         'ci','ci_lower','ci_upper',
                         'c_stat','c_stat_lb','c_stat_ub',
                         's0','avg_beta',
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
explore_categorical <- function(time,status,variable,dataset,path=getwd(),variables_list,
                                threshold,calib_quantiles=10,censor.t=5,make_plot=TRUE,
                                all_pop=TRUE){
  i <- 1
  out <- list()
  for (var in unique(data[[variable]])){
    subset <- dataset[get(variable)==var]
    n_event <- nrow(subset[get(status)==1]); n_total <- nrow(subset)
    total_pt <- sum(subset[[time]]); event_ir <- n_event/total_pt
    mod <- coxph(Surv(subset[[time]],subset[[status]]) ~ subset[[risk_score]])
    km <- survfit(Surv(subset[[time]],subset[[status]]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    
    # Refit model
    new_formula <- formula(paste0('Surv(subset[,get(time)],subset[,get(status)]) ~ ',paste0(variables_list,collapse=' + ')))
    mod <- coxph(new_formula,data=subset)
    subset$lp <- predict(mod,type='lp')
    
    # Calculate predicted risk with recalibration
    avg_beta <- mean(subset$lp)
    res <- coxph(Surv(subset[,get(time)],subset[,get(status)]) ~ lp, data=subset)
    km <- survfit(res, data=data.frame(x1=mean(lp)),type="kaplan-meier")
    s0 <- summary(km, times=c(5))$surv
    subset[,pred_risk := (1-(s0)^exp(lp - (avg_beta)))]
    
    # Tally binary test results
    tp <- nrow(subset[c(pred_risk >= threshold & get(status)==1)])
    tn <- nrow(subset[c(pred_risk < threshold & get(status)==0)])
    test_pos <- nrow(subset[pred_risk >= threshold]); test_neg <- nrow(subset[pred_risk < threshold])
    dz_pos <- nrow(subset[get(status)==1]); dz_neg <- nrow(subset[get(status)==0])
    subset[,':='(calib_groups = classifier(risk=subset$pred_risk,ncuts=calib_quantiles),
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
      incidence <- survivor(data=subset,risk_data="calib_groups",event=event,time=time,breakpoint=5)
      obv <- incidence$est*100
      pred <- subset[,mean(pred_risk),by='calib_groups']; setorder(pred,calib_groups)
      y_lim <- x_lim <- (max(obv,pred$V1*100,na.rm=T) %/% 5)*5+5
      pdf(file=paste0(path,'recalib_',variable,'_',var,'.pdf'),height=3,width=3,pointsize=3)
      plot(pred$V1*100,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
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
    out[[i]] <- data.frame(matrix(ncol=32,nrow=0))
    out[[i]] <- as.numeric(c(paste0(var),n_event,n_total,total_pt,event_ir,
                             ci,as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                             s0,avg_beta,
                             as.numeric(res$coefficients[1]),as.numeric(res$coefficients[1]-1.96*summary(res)$coefficients[3]),as.numeric(res$coefficients[1]+1.96*summary(res)$coefficients[3]),
                             relative_err,cum_err,gnd[2],gnd[3],
                             ppv,npv,sens,spec))
    names(out[[i]]) <- c('stratum','n_event','n_total','total_pt','event_ir',
                         'ci','ci_lower','ci_upper',
                         'c_stat','c_stat_lb','c_stat_ub',
                         's0','avg_beta',
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
    subset <- dataset
    n_event <- nrow(subset[get(status)==1]); n_total <- nrow(subset)
    total_pt <- sum(subset[[time]]); event_ir <- n_event/total_pt
    mod <- coxph(Surv(subset[[time]],subset[[status]]) ~ subset[[risk_score]])
    km <- survfit(Surv(subset[[time]],subset[[status]]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    
    # Refit model
    new_formula <- formula(paste0('Surv(subset[,get(time)],subset[,get(status)]) ~ ',paste0(variables_list,collapse=' + ')))
    mod <- coxph(new_formula,data=subset)
    subset$lp <- predict(mod,type='lp')
    
    # Calculate predicted risk with recalibration
    avg_beta <- mean(subset$lp)
    res <- coxph(Surv(subset[,get(time)],subset[,get(status)]) ~ lp, data=subset)
    km <- survfit(res, data=data.frame(x1=mean(lp)),type="kaplan-meier")
    s0 <- summary(km, times=c(5))$surv
    subset[,pred_risk := (1-(s0)^exp(lp - (avg_beta)))]
    
    # Binary test results
    tp <- nrow(subset[c(pred_risk >= threshold & get(status)==1)])
    tn <- nrow(subset[c(pred_risk < threshold & get(status)==0)])
    test_pos <- nrow(subset[pred_risk >= threshold]); test_neg <- nrow(subset[pred_risk < threshold])
    dz_pos <- nrow(subset[get(status)==1]); dz_neg <- nrow(subset[get(status)==0])
    subset[,':='(calib_groups = classifier(risk=subset$pred_risk,ncuts=calib_quantiles),
                 censored = ifelse(get(status)==0,1,0))]
    
    # GND component
    gnd_obj <- GND.calib(pred=subset[[pred_risk]],tvar=subset[[time]],out=subset[[status]],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      incidence <- survivor(data=subset,risk_data="calib_groups",event=event,time=time,breakpoint=5)
      obv <- incidence$est*100
      pred <- subset[,mean(pred_risk),by='calib_groups']; setorder(pred,calib_groups)
      y_lim <- x_lim <- (max(obv,pred$V1*100,na.rm=T) %/% 5)*5+5
      pdf(file=paste0(path,'recalib_all.pdf'),height=3,width=3,pointsize=3)
      plot(pred$V1*100,obv,xlab='Predicted',ylab='Observed',xlim=c(0,x_lim),ylim=c(0,y_lim),pch=19)
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
    out[[index]] <- data.frame(matrix(ncol=32,nrow=0))
    out[[index]] <- c('All',as.numeric(c(n_event,n_total,total_pt,event_ir,
                                         ci,as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                                         s0,avg_beta,
                                         as.numeric(res$coefficients[1]),as.numeric(res$coefficients[1]-1.96*summary(res)$coefficients[3]),as.numeric(res$coefficients[1]+1.96*summary(res)$coefficients[3]),
                                         relative_err,cum_err,gnd[2],gnd[3],
                                         ppv,npv,sens,spec)))
    names(out[[index]]) <- c('age','n_event','n_total','total_pt','event_ir',
                             'ci','ci_lower','ci_upper',
                             'c_stat','c_stat_lb','c_stat_ub',
                             's0','avg_beta',
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

# Scope covariate space
charge_vars <- c('start_fu_age_5','race_binary','ht_cm_atStartFu_10','wt_kg_atStartFu_15','sbp_atStartFu_20',
                 'dbp_atStartFu_10','tobacco_fin_prevAtstartFu','htn_fin_prevAtstartFu','dm_fin_prevAtstartFu',
                 'heartFailure_fin_prevAtstartFu','mi_fin_prevAtstartFu')

# Run explore functions
output_age_recal <- explore_continuous_recal(time='af_5y_sal.t',status='af_5y_sal',min_age=45,max_age=95,age_step=5,
                          variables_list=charge_vars,
                          age_variable='start_fu_age',dataset=vs,path='/data/arrhythmia/skhurshid/heterogeneity/',
                          threshold=0.05,make_plot=TRUE,all_pop=TRUE)
write.csv(output_age_recal,file='/data/arrhythmia/skhurshid/heterogeneity/charge_output_age_recal.csv')

output_sex_recal <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='Gender',
                                  variables_list=charge_vars,dataset=vs,
                                  path='/data/arrhythmia/skhurshid/heterogeneity/',threshold=0.05,make_plot=TRUE,all_pop=TRUE)
write.csv(output_sex_recal,file='/data/arrhythmia/skhurshid/heterogeneity/charge_output_sex_recal.csv')

output_hf_recal <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='heartFailure_fin_prevAtstartFu',
                                  variables_list=charge_vars,dataset=vs,
                                  path='/data/arrhythmia/skhurshid/heterogeneity/',threshold=0.05,make_plot=TRUE,all_pop=TRUE)
write.csv(output_hf_recal,file='/data/arrhythmia/skhurshid/heterogeneity/charge_output_hf_recal.csv')

output_stroke_recal <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='cvaTia_fin_prevAtstartFu',
                                 variables_list=charge_vars,dataset=vs,
                                 path='/data/arrhythmia/skhurshid/heterogeneity/',threshold=0.05,make_plot=TRUE,all_pop=TRUE)
write.csv(output_stroke_recal,file='/data/arrhythmia/skhurshid/heterogeneity/charge_output_stroke_recal.csv')

output_race_recal <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='race_binary',
                                     variables_list=charge_vars,dataset=vs,
                                     path='/data/arrhythmia/skhurshid/heterogeneity/',threshold=0.05,make_plot=TRUE,all_pop=TRUE)
write.csv(output_race_recal,file='/data/arrhythmia/skhurshid/heterogeneity/charge_output_race_recal.csv')
