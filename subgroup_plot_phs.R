# Script to generate subgroup analyses/plots in EHR-AF validation set

# Load validation set
load('/data/arrhythmia/skhurshid/ehr_af/vs_032120.RData')

# Dependenices
library(data.table)
library(survival)
library(timeROC)

setDT(vs)

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

# Define subgroups
f <- vs[Gender=='Female']
m <- vs[Gender=='Male']
w <- vs[race_binary==1]
nw <- vs[race_binary==2]
stroke <- vs[cvaTia_fin_prevAtstartFu==1]
hf <- vs[heartFailure_fin_prevAtstartFu==1]
big_dt <- list(f,m,w,nw,stroke,hf,vs)
  
# Obtain c-stat for each subgroup using timeROC bootstrap
auc <- list(); n <- 1
for (i in big_dt){
boot_score <- auc_boot(status='af_5y_sal',time='af_5y_sal.t',response='score',
                       data=i,times=4.9999,runs=200)
roc_score <-timeROC(T=i[,get('af_5y_sal.t')],delta=i[,get('af_5y_sal')],marker=i[,get('score')],
                    cause=1,weighting="marginal",times=4.999,iid=FALSE)$AUC[2]
auc[[n]] <- c(roc_score,roc_score-1.96*sd(boot_score),roc_score+1.96*sd(boot_score))
n <- n+1
}
auc_dt <- do.call(rbind,auc)

# Plot
pdf('/data/arrhythmia/skhurshid/broad_ibm_afrisk/subgroup_phs.pdf',height=2,width=4.4,
    pointsize=3)

par(oma=c(1,1,1,1))
par(mar=c(4,7.5,1,1))
col <- c('#e31a1c','#ff7f00','#a6d96a','#02818a','#3288bd','#5e4fa2','black')
plot(x=auc_dt[,1],y=seq(7,1,-1),xlim=c(0.600,0.825),ylim=c(0.75,7),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=c(rep(19,6),18),col=col,cex=c(rep(1.4,6),3),
     bty='n')
axis(1,at=seq(0.600,0.825,0.025),cex=1.4,pos=0)
axis(2,at=1:7,labels=FALSE,cex=1.4)
mtext('Subgroup',side=2,line=7,cex=1.3)
mtext('Concordance',side=1,line=3.5,cex=1.3)
segments(auc_dt[,2],7:1,auc_dt[,3],7:1,col=col,lwd=2)

plot_names <- c('Overall','Heart Failure','Stroke','Nonwhite','White','Male','Female')
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i]),
        side=2,line=1,las=2,at=i)
}

dev.off()
