library(ROCR); library(ggplot2); library(pROC); library(stringr); library(plyr); library(prodlim); library('rms'); library(data.table); library(survival);

DIRECTORY = "Y:/ibm-rwe_global/URI/Broad/output/"

FILE = 'Second_Data_Selection_Approach_AFib_Final_SUPERMART_122_Date_2019-10-04.csv'

setwd(DIRECTORY)
mydata <- read.csv(file = FILE, header = TRUE, sep = ',');

#mydata$num_days_index_date_to_AFib_combined_earliest_date[is.na(mydata$num_days_index_date_to_AFib_combined_earliest_date)] <- 1825

setDT(mydata)

mydata[,EHR_AF_Final_std := (EHR_AF_Final - (mean(EHR_AF_Final)))/sd((EHR_AF_Final))]
mydata[,CHARGE_AF_Final_std := (CHARGE_AF_Final - (mean(CHARGE_AF_Final)))/sd((CHARGE_AF_Final))]
mydata[,CHA2DS2_VASc_Final_std := (CHA2DS2_VASc_Final - (mean(CHA2DS2_VASc_Final)))/sd((CHA2DS2_VASc_Final))]
mydata[,C2HEST_Final_std := (C2HEST_Final - (mean(C2HEST_Final)))/sd((C2HEST_Final))]

###########################
gc()
memory.limit()
memory.limit(size = 16000)
memory.size()
###########################
## Estimated 5-yr AF risk from C2HEST study (PMID: 30292759)
chest_s0 <- exp(-0.005*5)
## Estimated C2HEST score from C2HEST study (PMID: 30292759)
chest_mean <- (88825*1 + 19270*2 + 8253 * 3 + 1373 * 4 + 90 * 5 + 15 * 6 + 15 * 7 + 15 * 8)/(45+90+1373+8253+19270+88825+310117)

## Estimated 5-yr AF risk from CHADSVASc study (PMID: 19762550)
chads_s0 <- exp(-0.023*5)
## Estimated C2HEST score from C2HEST study (PMID: 30292759)
chads_mean <- (162*1 + 184*2 + 203*3 + 208*4 + 95*5 + 57*6 + 25*7 + 9*8 + 1*9)/(1084)

### RISK SCORE
# Calculate pred 5 using published values
mydata[,pred_risk_uncalibrated_EHR_AF_Final := (1-0.9712209^exp(EHR_AF_Final - 6.728))*100]
mydata[,pred_risk_uncalibrated_CHARGE_AF_Final := (1-0.9718412736^exp(CHARGE_AF_Final - 12.58156))*100]
mydata[,pred_risk_uncalibrated_CHA2DS2_VASc_Final := (1-chads_s0^exp(CHA2DS2_VASc_Final - chads_mean))*100]
mydata[,pred_risk_uncalibrated_C2HEST_Final := (1-chest_s0^exp(C2HEST_Final - chest_mean))*100]

current_date = Sys.Date()

#######################################
### Calibration
#######################################
# Quantile sorter
classifier <- function(risk,ncuts){
  cuts <- quantile(risk,probs=seq(0,1,1/ncuts))
  index <- rep(NA,length(risk))
  for (i in 1:(length(cuts)-1)){
    for (j in 1:length(risk)){
      index[j] <- ifelse(risk[j] >= cuts[i],i,index[j])}}
  return(index)
}

#mydata$ehr_af_uncalibrated_dd <- classifier(risk=mydata$pred_risk_uncalibrated_EHR_AF_Final,ncuts=20)
mydata$charge_uncalibrated_dd <- classifier(risk=mydata$pred_risk_uncalibrated_CHARGE_AF_Final,ncuts=20)
#mydata$chest_uncalibrated <- round(mydata$pred_risk_uncalibrated_C2HEST_Final,2)
#mydata$chads_uncalibrated <- round(mydata$pred_risk_uncalibrated_CHA2DS2_VASc_Final,2)

# Function to generate survival estimates per AF risk quantile
survivor <- function(data,risk_data,time,status,eval.t){
  est <- rep(NA,times=length(unique(data[,risk_data])))
  lower <- rep(NA,times=length(unique(data[,risk_data])))
  upper <- rep(NA,times=length(unique(data[,risk_data])))
  for (i in 1:length(unique(data[,risk_data]))){
    km <- survfit(Surv(data[data[,risk_data]==unique(data[,risk_data])[order(unique(data[,risk_data]))][i],time],
                       data[data[,risk_data]==unique(data[,risk_data])[order(unique(data[,risk_data]))][i],status]) ~ 1)
    est[i] <- 1-stepfun(km$time, c(1, km$surv))(eval.t)
    upper[i] <- 1-stepfun(km$time, c(1, km$lower))(eval.t)
    lower[i] <- 1-stepfun(km$time, c(1, km$upper))(eval.t)
  }
  return(data.frame(est=est,upper=upper,lower=lower))
}

setDF(mydata)
#ehr_af_obv <- survivor(data=mydata,risk_data="ehr_af_uncalibrated_dd",time='num_days_to_AFib_or_to_cencor',status='is_AFib_combined',eval.t=365.25*5)
charge_obv <- survivor(data=mydata,risk_data="charge_uncalibrated_dd",time='num_days_to_AFib_or_to_cencor',status='is_AFib_combined',eval.t=365.25*5)
#chest_obv <- survivor(data=mydata,risk_data="chest_uncalibrated",time='num_days_to_AFib_or_to_cencor',status='is_AFib_combined',eval.t=365.25*5)
#chads_obv <- survivor(data=mydata,risk_data="chads_uncalibrated",time='num_days_to_AFib_or_to_cencor',status='is_AFib_combined',eval.t=365.25*5)
setDT(mydata)

### CALCULATE AVERAGE PREDICTED RISK IN EACH QUANTILE
#ehr_af_pred <- mydata[,mean(pred_risk_uncalibrated_EHR_AF_Final),by="ehr_af_uncalibrated_dd"][order(ehr_af_uncalibrated_dd)]
charge_pred <- mydata[,mean(pred_risk_uncalibrated_CHARGE_AF_Final),by="charge_uncalibrated_dd"][order(charge_uncalibrated_dd)]
#chest_pred <- mydata[,mean(pred_risk_uncalibrated_CHARGE_AF_Final),by="chest_uncalibrated"][order(chest_uncalibrated)]
#chads_pred <- mydata[,mean(pred_risk_uncalibrated_CHARGE_AF_Final),by="chads_uncalibrated"][order(chads_uncalibrated)]

# Generate the plots
## EHR-AF
pdf(file='~/cal_qual_charge.pdf',height=4,width=5,
    pointsize=3)
par(oma=c(2,3,1,1))
par(oma=c(2,3,1,1))

x <- charge_pred$V1
y <- charge_obv$est*100

plot(x,y,xlab='',ylab='',yaxt='n',
     xaxt='n',xlim=c(0,50),ylim=c(0,50),pch=19,cex=1.5)

axis(1,at=charge_pred$V1,
     labels=round(charge_pred$V1,1),cex.axis=1.7)
axis(2,cex.axis=1.6,at=seq(0,50,5),las=1)

segments(-1,-1,51,51,lwd=1.2,lty=2)

mtext("CHARGE-AF",side=3,cex=3,line=2,at=13.5)
mtext("5-yr predicted risk of AF (%)",side=1,cex=1.7,line=3)
mtext("Cumulative incidence of AF at 5 years (%)",side=2,cex=1.7,line=4.5)

dev.off()



