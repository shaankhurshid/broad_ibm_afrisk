# Dependencies
library(data.table); library(survival)

# Script to explore coefficients in subgroups

# Define explore function
explore_age <- function(time,status,age_variable,min_age,max_age,age_step,data){
  i <- 1
  out <- list()
  for (age in seq(min_age,max_age,age_step)){
    if(age == max_age-age_step){
      subset <- data[c(data[,age_variable] >= age & data[,age_variable] <= age+age_step),]
      n_af <- nrow(subset[subset[,status]==1,]); n_total <- nrow(subset)
      mod <- coxph(Surv(subset[,time],subset[,status]) ~ 
                     Gender + race_binary + age_10 + age_10_sq +
                     tobacco_fin_prevAtstartFu + ht_cm_atStartFu_10 + ht_10_sq +
                     wt_kg_atStartFu_15 + wt_15_sq + dbp80 + htn_fin_prevAtstartFu + 
                     lipid_fin_prevAtstartFu + heartFailure_fin_prevAtstartFu + cad_fin_prevAtstartFu +
                     valvDz_fin_prevAtstartFu + cvaTia_fin_prevAtstartFu + pad_fin_prevAtstartFu +
                     ckd_fin_prevAtstartFu + hypothyroid_fin_prevAtstartFu,data=subset)
      lp <- predict(mod, newdata=subset, type="lp")
      fit <- coxph(Surv(subset[,time],subset[,status])~lp,data=subset)
      out[[i]] <- data.frame(matrix(ncol=length(mod$coefficients)+9,nrow=0))
      out[[i]] <- rbind(out[[i]],c(paste0(age,'-',age+age_step-1),n_af,n_total,as.numeric(mod$coefficients),
                                   as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                                   as.numeric(fit$coefficients[1]),as.numeric(fit$coefficients[1]-1.96*summary(fit)$coefficients[3]),as.numeric(fit$coefficients[1]+1.96*summary(fit)$coefficients[3])),
                        c(paste0(age,'-',age+age_step-1),n_af,n_total,as.numeric(exp(mod$coefficients)),
                          as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                          as.numeric(fit$coefficients[1]),as.numeric(fit$coefficients[1]-1.96*summary(fit)$coefficients[3]),as.numeric(fit$coefficients[1]+1.96*summary(fit)$coefficients[3])))
      names(out[[i]]) <- c('age_class','n_af','n_total',names(mod$coefficients),'c_stat','c_stat_lb','c_stat_ub','cal','cal_lb','cal_ub')
      print(paste0('Just finished model ',i,' out of ',((max_age-min_age)/age_step),'!'))
      break
    } else {
      subset <- data[c(data[,age_variable] >= age & data[,age_variable] < age+age_step),]
      n_af <- nrow(subset[subset[,status]==1,]); n_total <- nrow(subset)
      mod <- coxph(Surv(subset[,time],subset[,status]) ~ 
                     Gender + race_binary + age_10 + age_10_sq +
                     tobacco_fin_prevAtstartFu + ht_cm_atStartFu_10 + ht_10_sq +
                     wt_kg_atStartFu_15 + wt_15_sq + dbp80 + htn_fin_prevAtstartFu + 
                     lipid_fin_prevAtstartFu + heartFailure_fin_prevAtstartFu + cad_fin_prevAtstartFu +
                     valvDz_fin_prevAtstartFu + cvaTia_fin_prevAtstartFu + pad_fin_prevAtstartFu +
                     ckd_fin_prevAtstartFu + hypothyroid_fin_prevAtstartFu,data=subset)
      lp <- predict(mod, newdata=subset, type="lp")
      fit <- coxph(Surv(subset[,time],subset[,status])~lp,data=subset)
      out[[i]] <- data.frame(matrix(ncol=length(mod$coefficients)+7,nrow=0))
      out[[i]] <- rbind(out[[i]],c(paste0(age,'-',age+age_step-1),n_af,n_total,as.numeric(mod$coefficients),
                                   as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                                   as.numeric(fit$coefficients[1]),as.numeric(fit$coefficients[1]-1.96*summary(fit)$coefficients[3]),as.numeric(fit$coefficients[1]+1.96*summary(fit)$coefficients[3])),
                        c(paste0(age,'-',age+age_step-1),n_af,n_total,as.numeric(exp(mod$coefficients)),
                          as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                          as.numeric(fit$coefficients[1]),as.numeric(fit$coefficients[1]-1.96*summary(fit)$coefficients[3]),as.numeric(fit$coefficients[1]+1.96*summary(fit)$coefficients[3])))
      names(out[[i]]) <- c('age_class','n_af','n_total',names(mod$coefficients),'c_stat','c_stat_lb','c_stat_ub','cal','cal_lb','cal_ub')
      print(paste0('Just finished model ',i,' out of ',((max_age-min_age)/age_step),'!'))
      i <- i+1
    }
  }
  return(do.call(rbind,out))
}

# Run explore function
data <- vs # put your time variable in time, AF variable in status, baseline age variable in age_variable

output <- explore_age(time='af_5y_sal.t',status='af_5y_sal',min_age=45,max_age=95,age_step=5,
                      age_variable='start_fu_age',data=data)

# Fix numerics as factors
setDT(output)
names <- names(output)[-1]
for (j in names){set(output,j=j,value=as.numeric(paste(output[[j]])))}
setDF(output)

# Output coefficients of main model
all_model <- coxph(Surv(data[,'af_5y_sal.t'],data[,'af_5y_sal']) ~ 
                     Gender + race_binary + age_10 + age_10_sq +
                     tobacco_fin_prevAtstartFu + ht_cm_atStartFu_10 + ht_10_sq +
                     wt_kg_atStartFu_15 + wt_15_sq + dbp80 + htn_fin_prevAtstartFu + 
                     lipid_fin_prevAtstartFu + heartFailure_fin_prevAtstartFu + cad_fin_prevAtstartFu +
                     valvDz_fin_prevAtstartFu + cvaTia_fin_prevAtstartFu + pad_fin_prevAtstartFu +
                     ckd_fin_prevAtstartFu + hypothyroid_fin_prevAtstartFu,data=data)
lp <- predict(all_model, newdata=data, type="lp")
fit <- coxph(Surv(data[,'af_5y_sal.t'],data[,'af_5y_sal'])~lp,data=data)
overall_coefs <- data.frame(matrix(ncol=length(all_model$coefficients)+7,nrow=0))
overall_coefs <- rbind(overall_coefs,c('all',sum(data$af_5y_sal),nrow(data),as.numeric(all_model$coefficients),
                                       as.numeric(summary(all_model)$concordance[1]),as.numeric(summary(all_model)$concordance[1]-1.96*summary(all_model)$concordance[2]),as.numeric(summary(all_model)$concordance[1]+1.96*summary(all_model)$concordance[2]),
                                       as.numeric(fit$coefficients[1]),as.numeric(fit$coefficients[1]-1.96*summary(fit)$coefficients[3]),as.numeric(fit$coefficients[1]+1.96*summary(fit)$coefficients[3])),
                       c('all',sum(data$af_5y_sal),nrow(data),exp(as.numeric(all_model$coefficients)),
                         as.numeric(summary(all_model)$concordance[1]),as.numeric(summary(all_model)$concordance[1]-1.96*summary(all_model)$concordance[2]),as.numeric(summary(all_model)$concordance[1]+1.96*summary(all_model)$concordance[2]),
                         as.numeric(fit$coefficients[1]),as.numeric(fit$coefficients[1]-1.96*summary(fit)$coefficients[3]),as.numeric(fit$coefficients[1]+1.96*summary(fit)$coefficients[3])))
names(overall_coefs) <- c('age_class','n_af','n_total',names(all_model$coefficients),'c_stat','c_stat_lb','c_stat_ub','cal','cal_lb','cal_ub')

# Fix numerics as factors
setDT(overall_coefs)
names <- names(overall_coefs)[-1]
for (j in names){set(overall_coefs,j=j,value=as.numeric(paste(overall_coefs[[j]])))}
setDF(overall_coefs)

# Combine overall results with subgroups results
output <- rbind(output,overall_coefs)

write.csv(output,file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/explore_beta_age.csv')

##################### Step 2: Graph coefficients for each model

# Graph
for (i in names(output)[4:22]){
  pdf(file=paste0('/data/arrhythmia/skhurshid/broad_ibm_afrisk/hr_plots/',i,'ageinc.pdf'),pointsize=3,
      height=3,width=3.5)
  plot(x=1:(length(output$age_class)/2),y=as.numeric(output[c(22,seq(2,20,2)),i]),pch=19,
       xaxt='n',yaxt='n',bty='n',xlab='',ylab='',col=c('red',rep('black',10)))
  axis(1,at=1:(length(output$age_class)/2),labels=unique(output$age_class)[c(11,1:10)])
  axis(2)
  mtext(1,line=3,text='age class')
  mtext(2,line=3,text='hazard ratio')
  mtext(3,line=0,text=paste(i))
  dev.off()
}
