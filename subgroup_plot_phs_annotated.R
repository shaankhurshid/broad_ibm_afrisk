# Script to generate subgroup analyses/plots in EHR-AF validation set

# Load validation set
load('/data/arrhythmia/skhurshid/ehr_af/vs_021219.RData')

# Dependenices
library(data.table)
library(survival)

setDT(vs)

# Define subgroups
f <- vs[Gender=='Female']
m <- vs[Gender=='Male']
hf <- vs[heartFailure_fin_prevAtstartFu==1]
stroke <- vs[cvaTia_fin_prevAtstartFu==1]

# Create all
c_f <- c(summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=f))$concordance[1],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=f))$concordance[1] - 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=f))$concordance[1] + 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         nrow(f),sum(f$af_5y_sal))
c_m <- c(summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=m))$concordance[1],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=m))$concordance[1] - 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=m))$concordance[1] + 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         nrow(m),sum(m$af_5y_sal))
c_stroke <- c(summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=stroke))$concordance[1],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=stroke))$concordance[1] - 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=stroke))$concordance[1] + 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         nrow(stroke),sum(stroke$af_5y_sal))
c_hf <- c(summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=hf))$concordance[1],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=hf))$concordance[1] - 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=hf))$concordance[1] + 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         nrow(hf),sum(hf$af_5y_sal))
c_overall <- c(summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[1],
          summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[1] - 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
          summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[1] + 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
          nrow(vs),sum(vs$af_5y_sal))

# bind
c_all <- rbind(c_f,c_m,c_stroke,c_hf,c_overall)

# Plot
pdf('/data/arrhythmia/skhurshid/broad_ibm_afrisk/subgroup_phs_annotated.pdf',height=1.5,width=4.4,
    pointsize=3)

par(oma=c(1,1,1,1))
par(mar=c(4,16.5,1,1))
col <- c('#e31a1c','#ff7f00','#33a02c','#1f78b4','black')
plot(x=c_all[,1],y=seq(5,1,-1),xlim=c(0.600,0.805),ylim=c(0.75,5.25),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=c(rep(19,4),18),col=col,cex=c(rep(1.4,4),3),
     bty='n')
axis(1,at=seq(0.600,0.825,0.025),cex=1.4,pos=0)
axis(2,at=1:18,labels=FALSE,cex=1.4)
mtext('Subgroup',side=2,line=15.5,cex=1.3)
mtext('Concordance',side=1,line=3.5,cex=1.3)
segments(c_all[,2],5:1,c_all[,3],5:1,col=col,lwd=2)

plot_names <- c('Overall','Heart Failure','Stroke','Male','Female')
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i],' (N=',c_all[6-i,4],'; AF=',c_all[6-i,5],')'),
        side=2,line=1,las=2,at=i)}
dev.off()
