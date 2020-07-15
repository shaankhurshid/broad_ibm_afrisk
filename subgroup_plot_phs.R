# Script to generate subgroup analyses/plots in EHR-AF validation set

# Load validation set
load('/data/arrhythmia/skhurshid/ehr_af/vs_032120.RData')

# Dependenices
library(data.table)
library(survival)

setDT(vs)

# Define subgroups
f <- vs[Gender=='Female']
m <- vs[Gender=='Male']
w <- vs[race_binary==1]
nw <- vs[race_binary==2]
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
c_w <- c(summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=w))$concordance[1],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=w))$concordance[1] - 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=w))$concordance[1] + 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         nrow(w),sum(w$af_5y_sal))
c_nw <- c(summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=nw))$concordance[1],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=nw))$concordance[1] - 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=nw))$concordance[1] + 1.96*summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=vs))$concordance[2],
         nrow(nw),sum(nw$af_5y_sal))
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
c_all <- rbind(c_f,c_m,c_w,c_nw,c_stroke,c_hf,c_overall)

# Plot
pdf('/data/arrhythmia/skhurshid/broad_ibm_afrisk/subgroup_phs.pdf',height=2,width=4.4,
    pointsize=3)

par(oma=c(1,1,1,1))
par(mar=c(4,7.5,1,1))
col <- c('#e31a1c','#ff7f00','#a6d96a','#02818a','#3288bd','#5e4fa2','black')
plot(x=c_all[,1],y=seq(7,1,-1),xlim=c(0.600,0.825),ylim=c(0.75,7),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=c(rep(19,6),18),col=col,cex=c(rep(1.4,6),3),
     bty='n')
axis(1,at=seq(0.600,0.825,0.025),cex=1.4,pos=0)
axis(2,at=1:7,labels=FALSE,cex=1.4)
mtext('Subgroup',side=2,line=7,cex=1.3)
mtext('Concordance',side=1,line=3.5,cex=1.3)
segments(c_all[,2],7:1,c_all[,3],7:1,col=col,lwd=2)

plot_names <- c('Overall','Heart Failure','Stroke','Nonwhite','White','Male','Female')
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i]),
        side=2,line=1,las=2,at=i)
}

dev.off()
