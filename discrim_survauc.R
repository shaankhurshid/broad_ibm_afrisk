# Script for calculating c-indices using survAUC package

# Dependencies
library(data.table)
library(survAUC)

# Load dataset
load(file='/Volumes/arrhythmia/skhurshid/ehr_af/vs_032120.RData')
setDT(vs)

data = vs[sample(1:nrow(vs),size=20000)]

# AUC CD
surv_obj <- Surv(data$af_5y_sal.t,data$af_5y_sal)
auc_chads <- AUC.cd(Surv.rsp=surv_obj,lp=data$chads,lpnew=data$chads,times=5)
auc_charge <- AUC.cd(Surv.rsp=surv_obj,lp=data$charge,lpnew=data$charge,times=5)
auc_chest <- AUC.cd(Surv.rsp=surv_obj,lp=data$chest,lpnew=data$chest,times=5)
auc_ehraf <- AUC.cd(Surv.rsp=surv_obj,lp=data$score,lpnew=data$score,times=5)

# Bootstrap function
auc_boot <- function(status,time,data,response,times=5,runs=200){
  out <- list()
  for (i in 1:runs){
    replicate <- data[sample(1:nrow(data),size=nrow(data),replace=TRUE)]
    surv_obj <- Surv(replicate[,get(time)],replicate[,get(status)])
    out[[i]] <- AUC.cd(Surv.rsp=surv_obj,lp=replicate[,get(response)],lpnew=replicate[,get(response)],times=times)
    if (i %% 50 == 0 ){print(paste0("I just finished ",i," out of ",runs," runs!"))}
  }
  return(do.call(rbind,out))
}

auc = auc_boot(status='af_5y_sal',time='af_5y_sal.t',response='score',data=a,runs=3)

# AUC CD
surv_obj <- Surv(vs$af_5y_sal.t,vs$af_5y_sal)
auc <- AUC.cd(Surv.rsp=surv_obj,lp=vs$score,lpnew=vs$score,times=5)

# Bootstrap function
auc_boot <- function(status,time,data,response,times=5,runs=200){
  out <- list()
  for (i in 1:runs){
  replicate <- data[sample(1:nrow(data),size=nrow(data),replace=TRUE)]
  surv_obj <- Surv(replicate[,get(time)],replicate[,get(status)])
  out[[i]] <- AUC.cd(Surv.rsp=surv_obj,lp=replicate[,get(response)],lpnew=replicate[,get(response)],times=times)
  if (i %% 50 == 0 ){print(paste0("I just finished ",i," out of ",runs," runs!"))}
  }
return(do.call(rbind,out))
}

auc = auc_boot(status='af_5y_sal',time='af_5y_sal.t',response='score',data=a,runs=3)