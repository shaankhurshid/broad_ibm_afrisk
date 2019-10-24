# Script to assess performance of EHR-AF in subgroups of interest

# Open latest validation set
load(file='/data/arrhythmia/data/ttac/shaan_rfiles/olivia/vs_021219.RData')

# Dependencies
library(rms)
library(data.table)

# Bootstrap function
boot <- function(time,status,response,data,runs,size){
  out <- rep(NA,times=runs)
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE),]
    out[i] <- concordance(coxph(Surv(sample[,time],sample[,status]) ~ sample[,response],data=sample))$concordance
    print(paste0('run ',i,' complete'))
  }
  return(out)
}

## AF
# EHR-AF ('time' = censored survival time, 'status' = outcome status, 'response' = classifier) by SEX
boot_score_f <- boot(time='af_5y_sal.t',status='af_5y_sal',response='score',data=vs[vs$Gender=='Female',],runs=200,size=nrow(vs[vs$Gender=='Female',]))
boot_score_m <- boot(time='af_5y_sal.t',status='af_5y_sal',response='score',data=vs[vs$Gender=='Male',],runs=200,size=nrow(vs[vs$Gender=='Male',]))

# Final AUC and 95% CI (1st value = AUC, 2nd value = Lower bound of 95% CI, 3rd value = Upper bound of 95% CI)
auc_score_f <- c(mean(boot_score_f),mean(boot_score_f)-1.96*sd(boot_score_f),mean(boot_score_f)+1.96*sd(boot_score_f))
auc_score_m <- c(mean(boot_score_m),mean(boot_score_m)-1.96*sd(boot_score_m),mean(boot_score_m)+1.96*sd(boot_score_m))