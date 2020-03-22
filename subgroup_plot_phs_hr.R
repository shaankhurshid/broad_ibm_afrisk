# Script to generate subgroup analyses/plots showing hazard ratio in EHR-AF validation set
# With N and AF annotations

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

# Standardize score in subgroups and overall
f[,ehraf_std := (score-(mean(score)))/(sd(score))]
m[,ehraf_std := (score-(mean(score)))/(sd(score))]
hf[,ehraf_std := (score-(mean(score)))/(sd(score))]
stroke[,ehraf_std := (score-(mean(score)))/(sd(score))]
vs[,ehraf_std := (score-(mean(score)))/(sd(score))]

# Create all
c_f <- c(exp(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=f)$coefficients[1]),
         exp(confint(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=f),'ehraf_std')[1]),
         exp(confint(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=f),'ehraf_std')[2]),
         nrow(f),sum(f$af_5y_sal))
c_m <- c(exp(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=m)$coefficients[1]),
         exp(confint(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=m),'ehraf_std')[1]),
         exp(confint(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=m),'ehraf_std')[2]),
         nrow(m),sum(m$af_5y_sal))
c_stroke <- c(exp(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=stroke)$coefficients[1]),
         exp(confint(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=stroke),'ehraf_std')[1]),
         exp(confint(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=stroke),'ehraf_std')[2]),
         nrow(stroke),sum(stroke$af_5y_sal))
c_hf <- c(exp(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=hf)$coefficients[1]),
              exp(confint(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=hf),'ehraf_std')[1]),
              exp(confint(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=hf),'ehraf_std')[2]),
              nrow(hf),sum(hf$af_5y_sal))
c_overall <- c(exp(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=vs)$coefficients[1]),
          exp(confint(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=vs),'ehraf_std')[1]),
          exp(confint(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ ehraf_std,data=vs),'ehraf_std')[2]),
          nrow(vs),sum(vs$af_5y_sal))

# bind
c_all <- rbind(c_f,c_m,c_stroke,c_hf,c_overall)

# Plot
pdf('/data/arrhythmia/skhurshid/broad_ibm_afrisk/subgroup_phs_hr.pdf',height=1.5,width=4.4,
    pointsize=3)

par(oma=c(1,1,1,1))
par(mar=c(4,7.5,1,1))
col <- c('#e31a1c','#ff7f00','#33a02c','#1f78b4','black')
plot(x=c_all[,1],y=seq(5,1,-1),xlim=c(1,3.50),ylim=c(0.75,5),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=c(rep(19,4),18),col=col,cex=c(rep(1.4,4),3),
     bty='n')
axis(1,at=seq(1,3.50,0.25),cex=1.4,pos=0)
axis(2,at=1:5,labels=FALSE,cex=1.4)
mtext('Subgroup',side=2,line=7,cex=1.3)
mtext('Standardized Hazard Ratio',side=1,line=3.5,cex=1.3)
segments(c_all[,2],5:1,c_all[,3],5:1,col=col,lwd=2)

plot_names <- c('Overall','Heart Failure','Stroke','Male','Female')
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i]),
        side=2,line=1,las=2,at=i)}
dev.off()
