# Script to generate subgroup analyses/plots showing hazard ratio in EHR-AF validation set
# With N and AF annotations

# Load validation set
load('/data/arrhythmia/skhurshid/ehr_af/vs_032120.RData')

# Dependenices
library(data.table)
library(survival)
library(timeROC)

# Bootstrap function
auc_boot <- function(status,time,data,response,times=5,runs=200){
  out <- list()
  for (i in 1:runs){
    replicate <- data[sample(1:nrow(data),size=nrow(data),replace=TRUE)]
    out[[i]] <- timeROC(T=replicate[,get(time)],delta=replicate[,get(status)],marker=replicate[,get(response)],
                        cause=1,weighting='marginal',times=times,iid=FALSE)$AUC[2]
    if (i %% 50 == 0 ){print(paste0("I just finished ",i," out of ",runs," runs!"))}
  }
  return(do.call(rbind,out))
}

boot_score <- auc_boot(status='af_5y_sal',time='af_5y_sal.t',response='score',
                       data=vs,times=4.9999,runs=3)
roc_score <-timeROC(T=vs$af_5y_sal.t,delta=vs$af_5y_sal,marker=vs$score,
                    cause=1,weighting="marginal",times=4.999,iid=FALSE)$AUC[2]
auc_score <- c(roc_score,roc_score-1.96*sd(boot_score),roc_score+1.96*sd(boot_score))

boot_chargeaf <- auc_boot(status='af_5y_sal',time='af_5y_sal.t',response='chargeaf',
                          data=vs,times=4.9999,runs=3)
roc_chargeaf <-timeROC(T=vs$af_5y_sal.t,delta=vs$af_5y_sal,marker=vs$chargeaf,
                       cause=1,weighting="marginal",times=4.999,iid=FALSE)$AUC[2]
auc_chargeaf <- c(roc_chargeaf,roc_chargeaf-1.96*sd(boot_chargeaf),roc_chargeaf+1.96*sd(boot_chargeaf))

boot_chest <- auc_boot(status='af_5y_sal',time='af_5y_sal.t',response='chest',
                       data=vs,times=4.9999,runs=3)
roc_chest <-timeROC(T=vs$af_5y_sal.t,delta=vs$af_5y_sal,marker=vs$chest,
                    cause=1,weighting="marginal",times=4.999,iid=FALSE)$AUC[2]
auc_chest <- c(roc_chest,roc_chest-1.96*sd(boot_chest),roc_chest+1.96*sd(boot_chest))

boot_chads <- auc_boot(status='af_5y_sal',time='af_5y_sal.t',response='chads',
                       data=vs,times=4.9999,runs=3)
roc_chads <-timeROC(T=vs$af_5y_sal.t,delta=vs$af_5y_sal,marker=vs$chads,
                    cause=1,weighting="marginal",times=4.999,iid=FALSE)$AUC[2]
auc_chads <- c(roc_chads,roc_chads-1.96*sd(boot_chads),roc_chads+1.96*sd(boot_chads))
