# Dependencies
library(data.table)
library(plyr)
library(survival)
library(nricens)

# Load dataset
load(file='/data/arrhythmia/skhurshid/ehr_af/vs_032120.RData')
load(file='/data/arrhythmia/skhurshid/ehr_af/ds_021219.RData')
setDT(vs); setDT(ds)

# Step 1: Create the CHARGE-AF OLD score in DS
## Model
charge_old <- coxph(Surv(af_5y_sal.t,af_5y_sal) ~ start_fu_age_5 + race_binary + ht_cm_atStartFu_10
                    + wt_kg_atStartFu_15 + sbp_atStartFu_20 + dbp_atStartFu_10
                    + tobacco_fin_prevAtstartFu + htn_fin_prevAtstartFu
                    + dm_fin_prevAtstartFu + heartFailure_fin_prevAtstartFu + mi_fin_prevAtstartFu,data=ds[start_fu_age >= 65])
## S0
# Mean beta in old set
charge_old_avgbeta <- mean(predict(charge_old,type='lp'))
# Linear predictor
ds[start_fu_age >= 65,old_charge_score := predict(charge_old,type='lp')]
# fit survival model
res <- coxph(Surv(af_5y_sal.t,af_5y_sal) ~ old_charge_score, data=ds[start_fu_age >= 65])
# set km as a survival function using the average level of each factor
km <- survfit(res, data=data.frame(x1=mean(old_charge_score)),type="kaplan-meier")
# set s0 as the survival coefficient of km
charge_old_s0 <- summary(km, times=c(5))$surv

# Step 2: Deploy among people on the cusp in the validation set
close <- vs[c(charge.pred5 >= 0.03 & charge.pred5 < 0.05 & start_fu_age >= 65)]
close[,old_charge_score := predict(charge_old,newdata=close,type='lp')]

### Calculate predicted risk
# Calculate pred 5 using survival equation
close[,charge_old_pred5 := (1-(charge_old_s0)^exp(old_charge_score - charge_old_avgbeta))*100]

### Step 3: Compare scores
# Quantile sorter
classifier <- function(risk,ncuts){
  cuts <- quantile(risk,probs=seq(0,1,1/ncuts))
  index <- rep(NA,length(risk))
  for (i in 1:(length(cuts)-1)){
    for (j in 1:length(risk)){
      index[j] <- ifelse(risk[j] >= cuts[i],i,index[j])}}
  return(index)
}

close$charge_decile <- classifier(risk=close$charge.pred5,ncuts=10)
close$charge_old_decile <- classifier(risk=close$charge_old_pred5,ncuts=10)

close$score_decile <- classifier(risk=close$pred5,ncuts=10)

# KM estimates per decile
survivor <- function(data,risk_data,event,time,breakpoint){
  af_est <- rep(NA,times=length(unique(data[,risk_data])))
  af_lower <- rep(NA,times=length(unique(data[,risk_data])))
  af_upper <- rep(NA,times=length(unique(data[,risk_data])))
  for (i in 1:length(unique(data[,risk_data]))){
    subset <- data[data[,risk_data]==unique(data[,risk_data])[order(unique(data[,risk_data]))][i],]
    km <- survfit(Surv(subset[,time],subset[,event]) ~ 1, data=subset)
    end_time <- which(round(km$time,1)==round(breakpoint,1))[1]
    af_est[i] <- 1-stepfun(km$time[1:end_time], c(1, km$surv[1:end_time]))(length(km$time[1:end_time]))
    af_upper[i] <- 1-stepfun(km$time[1:end_time], c(1, km$lower[1:end_time]))(length(km$time[1:end_time]))
    af_lower[i] <- 1-stepfun(km$time[1:end_time], c(1, km$upper[1:end_time]))(length(km$time[1:end_time]))
  }
  return(data.frame(af_est=af_est,af_upper=af_upper,af_lower=af_lower))
}

# Observed risks
setDF(close)
charge_obv <- survivor(data=close,risk_data="charge_decile",event='af_5y_sal',time='af_5y_sal.t',breakpoint=5.0)
charge_old_obv <- survivor(data=close,risk_data="charge_old_decile",event='af_5y_sal',time='af_5y_sal.t',breakpoint=5.0)
score_obv <- survivor(data=close,risk_data="score_decile",event='af_5y_sal',time='af_5y_sal.t',breakpoint=5.0)
setDT(close)

# Mean predicted risk per quantile
charge_pred <- close[,mean(charge.pred5)*100,by="charge_decile"][order(charge_decile)]
charge_old_pred <- close[,mean(charge_old_pred5),by="charge_old_decile"][order(charge_old_decile)]
score_pred <- close[,mean(pred5),by="score_decile"][order(score_decile)]

pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/close_charge_cal.pdf',height=5,width=5,
    pointsize=3)
par(oma=c(2,3,1,1))
par(oma=c(2,3,1,1))

x <- charge_pred$V1
y <- charge_obv[,1]*100

plot(x,y,xlab='',ylab='',yaxt='n',
     xaxt='n',xlim=c(0,10),ylim=c(0,10),pch=19,cex=1.5)

axis(1,at=seq(0,10,1),cex.axis=1.7)
axis(2,cex.axis=1.6,at=seq(0,10,1),las=1)

segments(-1,-1,11,11,lwd=1.2,lty=2)

mtext("CHARGE-AF Original",side=3,cex=1.7,line=2)
mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.7,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.7,line=4.5)

dev.off()

pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/close_charge_old_cal.pdf',height=5,width=5,
    pointsize=3)
par(oma=c(2,3,1,1))
par(oma=c(2,3,1,1))

x <- charge_old_pred$V1
y <- charge_old_obv[,1]*100

plot(x,y,xlab='',ylab='',yaxt='n',
     xaxt='n',xlim=c(0,10),ylim=c(0,10),pch=19,cex=1.5)

axis(1,at=seq(0,10,1),cex.axis=1.7)
axis(2,cex.axis=1.6,at=seq(0,10,1),las=1)

segments(-1,-1,51,51,lwd=1.2,lty=2)

mtext("CHARGE-AF Old",side=3,cex=1.7,line=2)
mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.7,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.7,line=4.5)

dev.off()

# Binary variables for high risk
close[,':='(charge_pred5_above5 = ifelse(charge.pred5 >= 0.05,1,0),
             charge_old_pred5_above5 = ifelse(charge_old_pred5 >= 5,1,0))]

charge_compare <- nricens(p.std=close$charge_pred5_above5, p.new=close$charge_old_pred5_above5,
                              time=close$af_5y_sal.t, event=close$af_5y_sal, cut=0.5,
                              niter = 10, t0=5,updown='category')