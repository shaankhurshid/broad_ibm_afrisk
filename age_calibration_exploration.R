# Dependencies
library(data.table)
library(plyr)
library(survival)
library(nricens)

# Load dataset
load(file='/data/arrhythmia/skhurshid/ehr_af/vs_032120.RData')
setDT(vs)

# Step 1: Find the age at which everyone has a predicted risk > 5% using CHARGE-AF
## First, find the CHARGE-AF score that = 5% risk in this dataset
charge_5 <- mean(vs[round(charge.pred5,4)==0.05]$chargeaf)
## Now, figure out what age is needed to get there assuming optimal risk factors otherwise
charge_optimal_noage <- (1*0.465 + mean(vs$ht_cm_atStartFu_10)*0.248 
                             + mean(vs$wt_kg_atStartFu_15)*0.115 + 6*0.197 + 8*(-0.101))
threshold_age <- (charge_5 - charge_optimal_noage)/0.1016 # 74.55
## Get everyone above that age
vs[,charge_age_threshold := ifelse(start_fu_age>=(threshold_age),1,0)]

# Assess metrics within this age
# Define explore function
explore_categorical <- function(time,status,variable,data,risk_score,risk,threshold){
  i <- 1
  out <- list()
  for (var in unique(data[,variable])){
    subset <- data[data[,variable]==unique(data[,variable])[i],]
    n_af <- nrow(subset[subset[,status]==1,]); n_total <- nrow(subset)
    mod <- coxph(Surv(subset[,time],subset[,status]) ~ subset[,risk_score],data=subset)
    km <- survfit(Surv(subset[,time],subset[,status]) ~ 1,data=subset)
    ci <- c((1-km$surv[length(km$surv)])*100,
            (1-km$upper[length(km$upper)])*100,
            (1-km$lower[length(km$lower)])*100)
    tp <- nrow(subset[c((subset[,risk] >= threshold) & (subset[,status]==1)),])
    tn <- nrow(subset[c((subset[,risk] < threshold) & (subset[,status]==0)),])
    test_pos <- nrow(subset[subset[,risk] >= threshold,]); test_neg <- nrow(subset[subset[,risk] < threshold,])
    dz_pos <- nrow(subset[subset[,status]==1,]); dz_neg <- nrow(subset[subset[,status]==0,])
    if (tp == 0 | tn == 0 | test_pos == 0 | dz_pos == 0){
      ppv <- c(tp/test_pos,NA,NA); npv <- c(tn/test_neg,NA,NA); sens <- c(tp/dz_pos,NA,NA); spec <- c(tn/dz_neg,NA,NA)
    } else {
      ppv <- c(tp/test_pos,binom.test(tp,test_pos)$conf.int[1],binom.test(tp,test_pos)$conf.int[2])
      npv <- c(tn/test_neg,binom.test(tn,test_neg)$conf.int[1],binom.test(tn,test_neg)$conf.int[2])
      sens <- c(tp/dz_pos,binom.test(tp,dz_pos)$conf.int[1],binom.test(tp,dz_pos)$conf.int[2])
      spec <- c(tn/dz_neg,binom.test(tn,dz_neg)$conf.int[1],binom.test(tn,dz_neg)$conf.int[2])
    }
    out[[i]] <- data.frame(matrix(ncol=24,nrow=0))
    out[[i]] <- c(paste0(unique(data[,variable])[i]),as.numeric(c(n_af,n_total,ci,
                                                                  as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                                                                  as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3]),
                                                                  ppv,npv,sens,spec)))
    names(out[[i]]) <- c('category','n_af','n_total','ci','ci_lower','ci_upper',
                         'c_stat','c_stat_lb','c_stat_ub','cal','cal_lb','cal_ub',
                         'ppv','ppv_lower','ppv_upper',
                         'npv','npv_lower','npv_upper',
                         'sens','sens_lower','sens_upper',
                         'spec','spec_lower','spec_upper')
    print(paste0('Just finished model ',i,' out of ',length(unique(data[,variable])),'!'))
    i <- i+1
  }
  return(data.frame(do.call(rbind,out)))
}

setDF(vs)
output_age <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='charge_age_threshold',
                                  data=vs,risk_score='chargeaf',risk='charge.pred5',threshold=0.05)
setDT(vs)

# Now re-derive the score in the older folks
vs_old <- vs[charge_age_threshold==1]

charge_old <- coxph(Surv(af_5y_sal.t,af_5y_sal) ~ start_fu_age_5 + race_binary + ht_cm_atStartFu_10
                    + wt_kg_atStartFu_15 + sbp_atStartFu_20 + dbp_atStartFu_10
                    + tobacco_fin_prevAtstartFu + htn_fin_prevAtstartFu
                    + dm_fin_prevAtstartFu + heartFailure_fin_prevAtstartFu + mi_fin_prevAtstartFu,data=vs_old)
vs_old[,charge_old := predict(charge_old,newdata=vs_old,type='lp')]

### Calculate predicted risk
# Average beta for charge
vs_old[,charge_avgbeta := mean(charge_old)]

# Fit survival model
res <- coxph(Surv(af_5y_sal.t,af_5y_sal) ~ charge_old, data=vs_old)

### derive s0
# set km as a survival function using the average level of each factor
km <- survfit(res, data=data.frame(x1=mean(charge_old)),type="kaplan-meier")
# Set s0 as the survival coefficient of km
charge_s0 <- summary(km, times=c(5))$surv

# Calculate pred 5 using survival equation
vs_old[,charge_old_pred5 := (1-(charge_s0)^exp(charge_old - charge_avgbeta))*100]

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

vs_old$charge_decile <- classifier(risk=vs_old$charge.pred5,ncuts=10)
vs_old$charge_old_decile <- classifier(risk=vs_old$charge_old_pred5,ncuts=10)

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
setDF(vs_old)
charge_obv <- survivor(data=vs_old,risk_data="charge_decile",event='af_5y_sal',time='af_5y_sal.t',breakpoint=5.0)
charge_old_obv <- survivor(data=vs_old,risk_data="charge_old_decile",event='af_5y_sal',time='af_5y_sal.t',breakpoint=5.0)
setDT(vs_old)

# Mean predicted risk per quantile
charge_pred <- vs_old[,mean(charge.pred5)*100,by="charge_decile"][order(charge_decile)]
charge_old_pred <- vs_old[,mean(charge_old_pred5),by="charge_old_decile"][order(charge_old_decile)]

pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/charge_cal.pdf',height=5,width=5,
    pointsize=3)
par(oma=c(2,3,1,1))
par(oma=c(2,3,1,1))

x <- charge_pred$V1
y <- charge_obv[,1]*100

plot(x,y,xlab='',ylab='',yaxt='n',
     xaxt='n',xlim=c(0,50),ylim=c(0,50),pch=19,cex=1.5)

axis(1,at=seq(0,50,10),cex.axis=1.7)
axis(2,cex.axis=1.6,at=seq(0,50,10),las=1)

segments(-1,-1,51,51,lwd=1.2,lty=2)

mtext("CHARGE-AF Original",side=3,cex=1.7,line=2)
mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.7,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.7,line=4.5)

dev.off()

pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/charge_old_cal.pdf',height=5,width=5,
    pointsize=3)
par(oma=c(2,3,1,1))
par(oma=c(2,3,1,1))

x <- charge_old_pred$V1
y <- charge_old_obv[,1]*100

plot(x,y,xlab='',ylab='',yaxt='n',
     xaxt='n',xlim=c(0,50),ylim=c(0,50),pch=19,cex=1.5)

axis(1,at=seq(0,50,10),cex.axis=1.7)
axis(2,cex.axis=1.6,at=seq(0,50,10),las=1)

segments(-1,-1,51,51,lwd=1.2,lty=2)

mtext("CHARGE-AF Old",side=3,cex=1.7,line=2)
mtext("Predicted risk of AF at 5 years (%)",side=1,cex=1.7,line=3)
mtext("Incidence of AF at 5 years (%)",side=2,cex=1.7,line=4.5)

dev.off()

# Binary variables for high risk
vs_old[,':='(charge_pred5_above5 = ifelse(charge.pred5 >= 0.05,1,0),
             charge_old_pred5_above5 = ifelse(charge_old_pred5 >= 5,1,0))]

# Oldies only
charge_compare <- nricens(p.std=vs_old$charge_pred5_above5, p.new=vs_old$charge_old_pred5_above5,
                              time=vs_old$af_5y_sal.t, event=vs_old$af_5y_sal, cut=0.5,
                              niter = 10, t0=5,updown='category')

