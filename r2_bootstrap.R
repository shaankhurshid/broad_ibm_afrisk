# Script to obtain Nagelkerke R2 with 95% CI using bootstrapping

# Dependencies
library(rms)

# Bootstrap function
boot <- function(time,status,response,data,runs,size){
  out <- rep(NA,times=runs)
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE),]
    out <- cph(Surv(sample[,time],sample[,status]) ~ sample[,response],data=sample)$stats['R2']
    print(paste0('run ',i,' complete'))
  }
  return(out)
}

## AF
# SCORE ('time' = censored survival time, 'status' = outcome status, 'response' = classifier)
actual_r2 <- cph(Surv(time,status) ~ score,data)$stats['R2']
boot_r2 <- boot(time='time',status='status',response='score',data=data,runs=200,size=nrow(data))

# Final AUC and 95% CI (1st value = AUC, 2nd value = Lower bound of 95% CI, 3rd value = Upper bound of 95% CI)
r2_score <- c(actual_r2,actual_r2-1.96*sd(boot_r2),actual_r2+1.96*sd(boot_r2))
