# Script to explore coefficients in subgroups

data <- vs

explore_age <- function(time,status,age_variable,min_age,max_age,age_step,data){
i <- 1
out <- list()
for (age in seq(min_age,max_age,age_step)){
  if(age == max_age-age_step){
  subset <- data[c(data[,age_variable] >= age & data[,age_variable] <= age+age_step),]
  mod <- coxph(Surv(subset[,time],subset[,status]) ~ 
                 Gender + age_10 + age_10_sq + race_binary +
                 tobacco_fin_prevAtstartFu + ht_cm_atStartFu_10 + ht_10_sq +
                 wt_kg_atStartFu_15 + wt_15_sq + dbp80 + htn_fin_prevAtstartFu + 
                 lipid_fin_prevAtstartFu + heartFailure_fin_prevAtstartFu + cad_fin_prevAtstartFu +
                 valvDz_fin_prevAtstartFu + cvaTia_fin_prevAtstartFu + pad_fin_prevAtstartFu +
                 ckd_fin_prevAtstartFu + hypothyroid_fin_prevAtstartFu,data=subset)
  lp <- predict(mod, newdata=subset, type="lp")
  fit <- coxph(Surv(subset[,time],subset[,status])~lp,data=subset)
  out[[i]] <- data.frame(matrix(ncol=length(mod$coefficients)+7,nrow=0))
  out[[i]] <- rbind(out[[i]],c(paste0(age,'-',age+age_step-1),mod$coefficients,
                               summary(mod)$concordance[1],summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2],summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2],
                               fit$coefficients[1],fit$coefficients[1]-1.96*summary(fit)$coefficients[3],fit$coefficients[1]+1.96*summary(fit)$coefficients[3]),
                    c(paste0(age,'-',age+age_step-1),exp(mod$coefficients),
                                                         summary(mod)$concordance[1],summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2],summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2],
                                                         fit$coefficients[1],fit$coefficients[1]-1.96*summary(fit)$coefficients[3],fit$coefficients[1]+1.96*summary(fit)$coefficients[3]))
  names(out[[i]]) <- c('age_class',names(mod$coefficients),'c_stat','c_stat_lb','c_stat_ub','cal','cal_lb','cal_ub')
  print(paste0('Just finished model ',i,' out of ',((max_age-min_age)/age_step),'!'))
  break
  } else {
    subset <- data[c(data[,age_variable] >= age & data[,age_variable] < age+age_step),]
    mod <- coxph(Surv(subset[,time],subset[,status]) ~ 
                   Gender + age_10 + age_10_sq + race_binary +
                   tobacco_fin_prevAtstartFu + ht_cm_atStartFu_10 + ht_10_sq +
                   wt_kg_atStartFu_15 + wt_15_sq + dbp80 + htn_fin_prevAtstartFu + 
                   lipid_fin_prevAtstartFu + heartFailure_fin_prevAtstartFu + cad_fin_prevAtstartFu +
                   valvDz_fin_prevAtstartFu + cvaTia_fin_prevAtstartFu + pad_fin_prevAtstartFu +
                   ckd_fin_prevAtstartFu + hypothyroid_fin_prevAtstartFu,data=subset)
    lp <- predict(mod, newdata=subset, type="lp")
    fit <- coxph(Surv(subset[,time],subset[,status])~lp,data=subset)
    out[[i]] <- data.frame(matrix(ncol=length(mod$coefficients)+7,nrow=0))
    out[[i]] <- rbind(out[[i]],c(paste0(age,'-',age+age_step-1),mod$coefficients,
                                 summary(mod)$concordance[1],summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2],summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2],
                                 fit$coefficients[1],fit$coefficients[1]-1.96*summary(fit)$coefficients[3],fit$coefficients[1]+1.96*summary(fit)$coefficients[3]),
                      c(paste0(age,'-',age+age_step-1),exp(mod$coefficients),
                                                           summary(mod)$concordance[1],summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2],summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2],
                                                           fit$coefficients[1],fit$coefficients[1]-1.96*summary(fit)$coefficients[3],fit$coefficients[1]+1.96*summary(fit)$coefficients[3]))
    names(out[[i]]) <- c('age_class',names(mod$coefficients),'c_stat','c_stat_lb','c_stat_ub','cal','cal_lb','cal_ub')
    print(paste0('Just finished model ',i,' out of ',((max_age-min_age)/age_step),'!'))
    i <- i+1
  }
}
return(do.call(rbind,out))
}

output <- explore_age(time='af_5y_sal.t',status='af_5y_sal',min_age=45,max_age=95,age_step=5,
                      age_variable='start_fu_age',data=data)

write.csv(output,file='/Volumes/arrhythmia/skhurshid/ehr_af/explore_beta_age.csv')

##################### Step 2: Graph coefficients for each model
# Output coefficients of main model
overall_coefs <- coxph(Surv(af_5y_sal.t,af_5y_sal) ~ 
               Gender + age_10 + age_10_sq + race_binary +
               tobacco_fin_prevAtstartFu + ht_cm_atStartFu_10 + ht_10_sq +
               wt_kg_atStartFu_15 + wt_15_sq + dbp80 + htn_fin_prevAtstartFu + 
               lipid_fin_prevAtstartFu + heartFailure_fin_prevAtstartFu + cad_fin_prevAtstartFu +
               valvDz_fin_prevAtstartFu + cvaTia_fin_prevAtstartFu + pad_fin_prevAtstartFu +
               ckd_fin_prevAtstartFu + hypothyroid_fin_prevAtstartFu,data=vs)$coefficients
overall_coefs <- rbind(overall_coefs,exp(overall_coefs))

# Graph
# 
