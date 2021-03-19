# Dependencies
rm(list=ls());
library(survival); library(data.table)
memory.limit(size = 456000)

# Source the suite (DF version)
source('subgroup_suite_df.R')

## Define explore function for continuous variables (input must be data.table)
explore_continuous_recal <- function(time,status,age_variable,min_age,max_age,variables_list,
                                     age_step,dataset,path=getwd(),threshold,
                                     calib_quantiles=10,censor.t=5,
                                     make_plot=TRUE,all_pop=TRUE,str1_file){
  
  
  i <- 1
  out <- list()
  for (age in seq(min_age,max_age-age_step,age_step)){
    
    if(age == max_age-age_step){
      subset <- dataset[c((dataset[,age_variable] >= age) & (dataset[,age_variable] <= age+age_step)),] 
    } else { 
      subset <- dataset[c((dataset[,age_variable] >= age) & (dataset[,age_variable] < age+age_step)),]}
    
    n_event <- nrow(subset[subset[,status]==1,]); n_total <- nrow(subset)
    total_pt <- sum(subset[,time]); event_ir <- n_event/total_pt
    km <- survfit(Surv(subset[,time],subset[,status]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    
    # Refit model
    print(paste0(ci))
    new_formula <- formula(paste0('Surv(subset[,time],subset[,status]) ~ ',paste0(variables_list,collapse=' + ')))
    mod <- coxph(new_formula,data=subset)
    subset$lp <- predict(mod,type='lp')
    
    # Calculate predicted risk with recalibration
    print('2')
    avg_beta <- mean(subset$lp)
    res <- coxph(Surv(subset[,time],subset[,status]) ~ lp, data=subset)
    km <- survfit(res, data=data.frame(x1=mean(lp)),type="kaplan-meier")
    s0 <- summary(km, times=c(censor.t))$surv
    print(paste0(s0))
    subset$pred_risk <- (1-(s0)^exp(subset$lp - (avg_beta)))
    
    # Sens/spec metrics
    print('3')
    tp <- nrow(subset[c((subset$pred_risk >= threshold) & subset[,status]==1),])
    tn <- nrow(subset[c((subset$pred_risk < threshold) & subset[,status]==0),])
    test_pos <- nrow(subset[subset$pred_risk >= threshold,]); test_neg <- nrow(subset[subset$pred_risk < threshold,])
    dz_pos <- nrow(subset[subset[,status]==1,]); dz_neg <- nrow(subset[subset[,status]==0,])
    subset$calib_groups <- classifier(risk=subset$pred_risk,ncuts=calib_quantiles)
    subset$censored <- ifelse(subset[,status]==0,1,0)
    
    # GND component
    gnd_obj <- GND.calib(pred=subset$pred_risk,tvar=subset[,time],out=subset[,status],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data="calib_groups",event=status,time=time,breakpoint=censor.t)
      obv <- incidence$est*100
      pred <- unlist(lapply(split(subset$pred_risk,subset$calib_groups),mean,na.rm=TRUE))
      y_lim <- x_lim <- (max(obv,pred*100,na.rm=T) %/% 5)*5+5
      #pdf(file=paste0(path,'calib_',age,'.pdf'),height=3,width=3,pointsize=3)
      pdf(file=paste0(path,paste0(str1_file,'calib_',age,'.pdf')),height=3,width=3,pointsize=3)
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
    n_event <- nrow(subset[subset[,status]==1,]); n_total <- nrow(subset)
    total_pt <- sum(subset[,time]); event_ir <- n_event/total_pt
    km <- survfit(Surv(subset[,time],subset[,status]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    
    # Refit model
    new_formula <- formula(paste0('Surv(subset[,time],subset[,status]) ~ ',paste0(variables_list,collapse=' + ')))
    mod <- coxph(new_formula,data=subset)
    subset$lp <- predict(mod,type='lp')
    
    # Calculate predicted risk with recalibration
    avg_beta <- mean(subset$lp)
    res <- coxph(Surv(subset[,time],subset[,status]) ~ lp, data=subset)
    km <- survfit(res, data=data.frame(x1=mean(lp)),type="kaplan-meier")
    s0 <- summary(km, times=c(censor.t))$surv
    subset$pred_risk <- (1-(s0)^exp(subset$lp - (avg_beta)))
    
    # Sens/spec metrics
    tp <- nrow(subset[c((subset$pred_risk >= threshold) & subset[,status]==1),])
    tn <- nrow(subset[c((subset$pred_risk < threshold) & subset[,status]==0),])
    test_pos <- nrow(subset[subset$pred_risk >= threshold,]); test_neg <- nrow(subset[subset$pred_risk < threshold,])
    dz_pos <- nrow(subset[subset[,status]==1,]); dz_neg <- nrow(subset[subset[,status]==0,])
    subset$calib_groups <- classifier(risk=subset$pred_risk,ncuts=calib_quantiles)
    subset$censored <- ifelse(subset[,status]==0,1,0)
    
    # GND component
    gnd_obj <- GND.calib(pred=subset$pred_risk,tvar=subset[,time],out=subset[,status],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data='calib_groups',event=status,time=time,breakpoint=censor.t)
      obv <- incidence$est*100
      pred <- unlist(lapply(split(subset$pred_risk,subset$calib_groups),mean,na.rm=TRUE))
      y_lim <- x_lim <- (max(obv,pred*100,na.rm=T) %/% 5)*5+5
      #pdf(file=paste0(path,'calib_',age,'.pdf'),height=3,width=3,pointsize=3)
      pdf(file=paste0(path,paste0(str1_file,'calib_all','.pdf')),height=3,width=3,pointsize=3)
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
explore_categorical_refit <- function(time,status,variable,dataset,path=getwd(),variables_list,
                                      threshold,calib_quantiles=10,censor.t,make_plot=TRUE,
                                      all_pop=TRUE,str1_file){
  i <- 1
  out <- list()
  for (var in unique(dataset[,variable])){
    subset <- dataset[dataset[,variable]==var,]
    n_event <- nrow(subset[subset[,status]==1,]); n_total <- nrow(subset)
    total_pt <- sum(subset[,time]); event_ir <- n_event/total_pt
    mod <- coxph(Surv(subset[,time],subset[,status]) ~ subset[,risk_score])
    km <- survfit(Surv(subset[,time],subset[,status]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    
    # Refit model
    new_formula <- formula(paste0('Surv(subset[,time],subset[,status]) ~ ',paste0(variables_list,collapse=' + ')))
    mod <- coxph(new_formula,data=subset)
    subset$lp <- predict(mod,type='lp')
    
    # Calculate predicted risk with recalibration
    avg_beta <- mean(subset$lp)
    res <- coxph(Surv(subset[,time],subset[,status]) ~ lp, data=subset)
    km <- survfit(res, data=data.frame(x1=mean(lp)),type="kaplan-meier")
    s0 <- summary(km, times=c(censor.t))$surv
    subset$pred_risk <- (1-(s0)^exp(subset$lp - (avg_beta)))
    
    # Sens/spec metrics
    tp <- nrow(subset[c((subset$pred_risk >= threshold) & subset[,status]==1),])
    tn <- nrow(subset[c((subset$pred_risk < threshold) & subset[,status]==0),])
    test_pos <- nrow(subset[subset$pred_risk >= threshold,]); test_neg <- nrow(subset[subset$pred_risk < threshold,])
    dz_pos <- nrow(subset[subset[,status]==1,]); dz_neg <- nrow(subset[subset[,status]==0,])
    subset$calib_groups <- classifier(risk=subset$pred_risk,ncuts=calib_quantiles)
    subset$censored <- ifelse(subset[,status]==0,1,0)
    
    # GND component
    gnd_obj <- GND.calib(pred=subset$pred_risk,tvar=subset[,time],out=subset[,status],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data='calib_groups',event=status,time=time,breakpoint=censor.t)
      obv <- incidence$est*100
      pred <- unlist(lapply(split(subset$pred_risk,subset$calib_groups),mean,na.rm=TRUE))
      y_lim <- x_lim <- (max(obv,pred*100,na.rm=T) %/% 5)*5+5
      #pdf(file=paste0(path,'calib_',var,'.pdf'),height=3,width=3,pointsize=3)
      pdf(file=paste0(path,paste0(str1_file,'calib_',var,'.pdf')),height=3,width=3,pointsize=3)
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
    n_event <- nrow(subset[subset[,status]==1,]); n_total <- nrow(subset)
    total_pt <- sum(subset[,time]); event_ir <- n_event/total_pt
    km <- survfit(Surv(subset[,time],subset[,status]) ~ 1)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    
    # Refit model
    new_formula <- formula(paste0('Surv(subset[,time],subset[,status]) ~ ',paste0(variables_list,collapse=' + ')))
    mod <- coxph(new_formula,data=subset)
    subset$lp <- predict(mod,type='lp')
    
    # Calculate predicted risk with recalibration
    avg_beta <- mean(subset$lp)
    res <- coxph(Surv(subset[,time],subset[,status]) ~ lp, data=subset)
    km <- survfit(res, data=data.frame(x1=mean(lp)),type="kaplan-meier")
    s0 <- summary(km, times=c(censor.t))$surv
    subset$pred_risk <- (1-(s0)^exp(subset$lp - (avg_beta)))
    
    # Sens/spec metrics
    tp <- nrow(subset[c((subset$pred_risk >= threshold) & subset[,status]==1),])
    tn <- nrow(subset[c((subset$pred_risk < threshold) & subset[,status]==0),])
    test_pos <- nrow(subset[subset$pred_risk >= threshold,]); test_neg <- nrow(subset[subset$pred_risk < threshold,])
    dz_pos <- nrow(subset[subset[,status]==1,]); dz_neg <- nrow(subset[subset[,status]==0,])
    subset$calib_groups <- classifier(risk=subset$pred_risk,ncuts=calib_quantiles)
    subset$censored <- ifelse(subset[,status]==0,1,0)
    
    # GND component
    gnd_obj <- GND.calib(pred=subset$pred_risk,tvar=subset[,time],out=subset[,status],groups=subset$calib_groups,
                         cens.t=subset$censored,adm.cens = censor.t)
    gnd <- gnd_obj[[2]]
    relative_err <- mean(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc)/gnd_obj[[1]]$kmperc)
    cum_err <- sum(abs(gnd_obj[[1]]$expectedperc-gnd_obj[[1]]$kmperc))
    
    # Plotting component
    if (make_plot == TRUE){
      # risk_data = stratum; event = desired event; time = time to event (IN YEARS); breakpoint = time to evaluate (IN YEARS)
      incidence <- survivor(data=subset,risk_data='calib_groups',event=status,time=time,breakpoint=censor.t)
      obv <- incidence$est*100
      pred <- unlist(lapply(split(subset$pred_risk,subset$calib_groups),mean,na.rm=TRUE))
      y_lim <- x_lim <- (max(obv,pred*100,na.rm=T) %/% 5)*5+5
      #pdf(file=paste0(path,'calib_',var,'.pdf'),height=3,width=3,pointsize=3)
      pdf(file=paste0(path,paste0(str1_file,'calib_all','.pdf')),height=3,width=3,pointsize=3)
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

#########################
### AF & CHARGE model ###
#########################

#DIRECTORY = "Y:/ibm-rwe_global/URI/Broad/Validation_Tool/output/"; setwd(DIRECTORY)
#ILE = 'AFib_Data_Selection_Final_SUPERMART_122_Date_2021-01-05.csv'
#s <- read.csv(file = FILE, header = TRUE, sep = ',');
vs <- fread(file='/data/arrhythmia/skhurshid/heterogeneity/af_test.csv')

setDT(vs)
names(vs)[names(vs) == "is_AFib_combined"] <- "af_5y_sal"
names(vs)[names(vs) == "num_days_to_AFib_or_to_cencor"] <- "af_5y_sal.t"
names(vs)[names(vs) == "CHARGE_AF_Final"] <- "chargeaf"
names(vs)[names(vs) == "age"] <- "start_fu_age"
names(vs)[names(vs) == "Gender"] <- "is_male"

#Converting time to event from days to years.
vs$af_5y_sal.t <- vs$af_5y_sal.t / 365
vs$af_5y_sal.t <- as.numeric(vs$af_5y_sal.t)

#DIRECTORY = "Y:/ibm-rwe_global/URI/Broad/"; setwd(DIRECTORY)

#vs[,score_avgbeta_chargeaf := mean(chargeaf)]
#res_chargeaf <- coxph(Surv(af_5y_sal.t,af_5y_sal) ~ chargeaf, data = vs); summary(res_chargeaf)
#km_chargeaf <- survfit(res_chargeaf, data=data.frame(x1 = mean(chargeaf)),type="kaplan-meier")
#s0_chargeaf <- summary(km_chargeaf, times = c(5))$surv; s0_chargeaf
#vs[,charge.pred5 := (1-(s0_chargeaf)^exp(chargeaf - (score_avgbeta_chargeaf)))]; summary(vs$charge.pred5)

setDF(vs)
#setDT(vs)

# Scope covariate space
charge_vars <- c('start_fu_age','is_white','is_male') # Replace with all when run for real

#charge_vars <- c('start_fu_age_5','race_binary','ht_cm_atStartFu_10','wt_kg_atStartFu_15','sbp_atStartFu_20',
#                'dbp_atStartFu_10','tobacco_fin_prevAtstartFu','htn_fin_prevAtstartFu','dm_fin_prevAtstartFu',
#                'heartFailure_fin_prevAtstartFu','mi_fin_prevAtstartFu')

### REPLACE ABOVE WITH CHARGE COMPONENT COLUMNS
output_age <- explore_continuous_recal(time='af_5y_sal.t',status='af_5y_sal',min_age=45,max_age=65,age_step=20,
                          age_variable='start_fu_age',data=vs,path='~/',variables_list=charge_vars,
                          threshold=0.05,censor.t=5,make_plot=TRUE,all_pop=TRUE,str1_file='Inc_AF_age_')
write.csv(output_age,file='Inc_AF_charge_output_age.csv')

# Run explore functions
output_age <- explore_age(time='af_5y_sal.t',status='af_5y_sal',min_age=45,max_age=90,age_step=5,
                          age_variable='start_fu_age',data=vs,path=WORK_DIR,variables_list=charge_vars,
                          risk_score='chargeaf',pred_risk='charge.pred5',threshold=0.05,censor.t=5,make_plot=TRUE,all_pop=TRUE,str1_file='Inc_AF_age_')
write.csv(output_age,file='Inc_AF_charge_output_age.csv')

output_sex <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='is_male',
                                  data=vs,risk_score='chargeaf',pred_risk='charge.pred5',variables_list=charge_vars,
                                  path=WORK_DIR,threshold=0.05,censor.t=5,make_plot=TRUE,all_pop=TRUE,str1_file='Inc_AF_is_male_')
write.csv(output_sex,file='Inc_AF_charge_output_sex.csv')

output_hf <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='Is_HF_ICD',
                                 data=vs,risk_score='chargeaf',pred_risk='charge.pred5',variables_list=charge_vars,
                                 path=WORK_DIR,threshold=0.05,censor.t=5,make_plot=TRUE,all_pop=TRUE,str1_file='Inc_AF_Is_HF_ICD_')
write.csv(output_hf,file='Inc_AF_charge_output_hf.csv')

output_stroke <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='Is_Stroke_TIA_ICD',
                                     data=vs,risk_score='chargeaf',pred_risk='charge.pred5',variables_list=charge_vars,
                                     path=WORK_DIR,threshold=0.05,censor.t=5,make_plot=TRUE,all_pop=TRUE,str1_file='Inc_AF_Is_Stroke_TIA_ICD_')
write.csv(output_stroke,file='Inc_AF_charge_output_stroke.csv')

output_race <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='is_white',
                                   data=vs,risk_score='chargeaf',pred_risk='charge.pred5',variables_list=charge_vars,
                                   path=WORK_DIR,threshold=0.05,censor.t=5,make_plot=TRUE,all_pop=TRUE,str1_file='Inc_AF_is_white_')
write.csv(output_race,file='Inc_AF_charge_output_race.csv')

#########################
### ASCVD & PCE model ###
#########################

DIRECTORY = "Y:/ibm-rwe_global/URI/Broad/Validation_Tool/output/"; setwd(DIRECTORY)
FILE = 'ASCVD_Data_Selection_Final_SUPERMART_122_Date_2021-01-05.csv'
vs <- read.csv(file = FILE, header = TRUE, sep = ',');

vs$num_days_to_ASCVD_or_to_cencor <- vs$num_days_to_ASCVD_or_to_cencor / 365
vs$num_days_to_ASCVD_or_to_cencor <- as.numeric(vs$num_days_to_ASCVD_or_to_cencor)

vs$ASCVD_Score_Percent <- vs$ASCVD_Score_Percent / 100

DIRECTORY = "Y:/ibm-rwe_global/URI/Broad/"; setwd(DIRECTORY)

vs_temp <- vs
vs <- vs_temp

# Scope covariate space
ascvd_vars <- c('start_fu_age_5','race_binary','ht_cm_atStartFu_10','wt_kg_atStartFu_15','sbp_atStartFu_20',
                 'dbp_atStartFu_10','tobacco_fin_prevAtstartFu','htn_fin_prevAtstartFu','dm_fin_prevAtstartFu',
                 'heartFailure_fin_prevAtstartFu','mi_fin_prevAtstartFu')
### REPLACE ABOVE WITH ASCVD COMPONENT COLUMNS

str2_file = ""; str3_file = "";
#vs <- vs[which(vs$is_male == 1 & vs$is_white == 1),]; str2_file = 'TABLE_A_2013_Male_White_'; str3_file = '_TABLE_A_2013_Male_White'
#vs <- vs[which(vs$is_male == 1 & vs$is_white == 0),]; str2_file = 'TABLE_A_2013_Male_Nonwhite_'; str3_file = '_TABLE_A_2013_Male_Nonwhite'
#vs <- vs[which(vs$is_male == 0 & vs$is_white == 1),]; str2_file = 'TABLE_A_2013_Female_White_'; str3_file = '_TABLE_A_2013_Female_White'
vs <- vs[which(vs$is_male == 0 & vs$is_white == 0),]; str2_file = 'TABLE_A_2013_Female_Nonwhite_'; str3_file = '_TABLE_A_2013_Female_Nonwhite'

setDF(vs)
# Run explore functions
output_age <- explore_age(time='num_days_to_ASCVD_or_to_cencor',status='OUTCOME_is_ASCVD_in_follow_up',min_age=45,max_age=80,age_step=5,
                          age_variable='age',data=vs,path=WORK_DIR,variables_list=ascvd_vars,
                          risk_score='Individual_Sum_for_ASCVD_Score',pred_risk='ASCVD_Score_Percent',threshold=0.075,censor.t=10,make_plot=TRUE,all_pop=TRUE,str1_file=paste0('Inc_ASCVD_age_', str2_file))
write.csv(output_age,file=paste0('Inc_ASCVD_pce_output_age', str3_file, '.csv'))

output_sex <- explore_categorical(time='num_days_to_ASCVD_or_to_cencor',status='OUTCOME_is_ASCVD_in_follow_up',variable='is_male',
                                  data=vs,risk_score='Individual_Sum_for_ASCVD_Score',pred_risk='ASCVD_Score_Percent',variables_list=ascvd_vars,
                                  path=WORK_DIR,threshold=0.075,censor.t=10,make_plot=TRUE,all_pop=TRUE,str1_file=paste0('Inc_ASCVD_is_male_', str2_file))
write.csv(output_sex,file=paste0('Inc_ASCVD_pce_output_sex', str3_file, '.csv'))

output_hf <- explore_categorical(time='num_days_to_ASCVD_or_to_cencor',status='OUTCOME_is_ASCVD_in_follow_up',variable='Is_HF_ICD',
                                 data=vs,risk_score='Individual_Sum_for_ASCVD_Score',pred_risk='ASCVD_Score_Percent',variables_list=ascvd_vars,
                                 path=WORK_DIR,threshold=0.075,censor.t=10,make_plot=TRUE,all_pop=TRUE,str1_file=paste0('Inc_ASCVD_Is_HF_ICD_', str2_file))
write.csv(output_hf,file=paste0('Inc_ASCVD_pce_output_hf', str3_file, '.csv'))

output_race <- explore_categorical(time='num_days_to_ASCVD_or_to_cencor',status='OUTCOME_is_ASCVD_in_follow_up',variable='is_white',
                                   data=vs,risk_score='Individual_Sum_for_ASCVD_Score',pred_risk='ASCVD_Score_Percent',variables_list=ascvd_vars,
                                   path=WORK_DIR,threshold=0.075,censor.t=10,make_plot=TRUE,all_pop=TRUE,str1_file=paste0('Inc_ASCVD_is_white_', str2_file))
write.csv(output_race,file=paste0('Inc_ASCVD_pce_output_race', str3_file, '.csv'))


