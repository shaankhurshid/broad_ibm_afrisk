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

# Model for standardized scores
mod_EHR_AF_Final <- coxph(Surv(num_days_to_AFib_or_to_cencor,is_AFib_combined) ~ EHR_AF_Final_std, data=mydata); summary(mod_EHR_AF_Final); AIC(mod_EHR_AF_Final)
mod_CHARGE_AF_Final <- coxph(Surv(num_days_to_AFib_or_to_cencor,is_AFib_combined) ~ CHARGE_AF_Final_std, data=mydata); summary(mod_CHARGE_AF_Final); AIC(mod_CHARGE_AF_Final)
mod_CHA2DS2_VASc_Final <- coxph(Surv(num_days_to_AFib_or_to_cencor,is_AFib_combined) ~ CHA2DS2_VASc_Final_std, data=mydata); summary(mod_CHA2DS2_VASc_Final); AIC(mod_CHA2DS2_VASc_Final)
mod_C2HEST_Final <- coxph(Surv(num_days_to_AFib_or_to_cencor,is_AFib_combined) ~ C2HEST_Final_std, data=mydata); summary(mod_C2HEST_Final); AIC(mod_C2HEST_Final)

### RISK SCORE
# Average beta
mydata[,score_avgbeta_EHR_AF_Final := mean(EHR_AF_Final)]
mydata[,score_avgbeta_CHARGE_AF_Final := mean(CHARGE_AF_Final)]
mydata[,score_avgbeta_CHA2DS2_VASc_Final := mean(CHA2DS2_VASc_Final)]
mydata[,score_avgbeta_C2HEST_Final := mean(C2HEST_Final)]

# Fit survival model
res_EHR_AF_Final <- coxph(Surv(num_days_to_AFib_or_to_cencor,is_AFib_combined) ~ EHR_AF_Final, data = mydata); summary(res_EHR_AF_Final)
res_CHARGE_AF_Final <- coxph(Surv(num_days_to_AFib_or_to_cencor,is_AFib_combined) ~ CHARGE_AF_Final, data = mydata); summary(res_CHARGE_AF_Final)
res_CHA2DS2_VASc_Final <- coxph(Surv(num_days_to_AFib_or_to_cencor,is_AFib_combined) ~ CHA2DS2_VASc_Final, data = mydata); summary(res_CHA2DS2_VASc_Final)
res_C2HEST_Final <- coxph(Surv(num_days_to_AFib_or_to_cencor,is_AFib_combined) ~ C2HEST_Final, data = mydata); summary(res_C2HEST_Final)

### derive s0
# set km as a survival function using the average level of each factor
km_EHR_AF_Final <- survfit(res_EHR_AF_Final, data=data.frame(x1 = mean(EHR_AF_Final)),type="kaplan-meier")
km_CHARGE_AF_Final <- survfit(res_CHARGE_AF_Final, data=data.frame(x1 = mean(CHARGE_AF_Final)),type="kaplan-meier")
km_CHA2DS2_VASc_Final <- survfit(res_CHA2DS2_VASc_Final, data=data.frame(x1 = mean(CHA2DS2_VASc_Final)),type="kaplan-meier")
km_C2HEST_Final <- survfit(res_C2HEST_Final, data=data.frame(x1 = mean(C2HEST_Final)),type="kaplan-meier")

# Set s0 as the survival coefficient of km
s0_EHR_AF_Final <- summary(km_EHR_AF_Final, times = c(5 * 365))$surv; s0_EHR_AF_Final
s0_CHARGE_AF_Final <- summary(km_CHARGE_AF_Final, times = c(5 * 365))$surv; s0_CHARGE_AF_Final
s0_CHA2DS2_VASc_Final <- summary(km_CHA2DS2_VASc_Final, times = c(5 * 365))$surv; s0_CHA2DS2_VASc_Final
s0_C2HEST_Final <- summary(km_C2HEST_Final, times = c(5 * 365))$surv; s0_C2HEST_Final

# Calculate pred 5 using survival equation
mydata[,pred_risk_EHR_AF_Final := (1-(s0_EHR_AF_Final)^exp(EHR_AF_Final - (score_avgbeta_EHR_AF_Final)))*100]; summary(mydata$pred_risk_EHR_AF_Final)
mydata[,pred_risk_CHARGE_AF_Final := (1-(s0_CHARGE_AF_Final)^exp(CHARGE_AF_Final - (score_avgbeta_CHARGE_AF_Final)))*100]; summary(mydata$pred_risk_CHARGE_AF_Final)
mydata[,pred_risk_CHA2DS2_VASc_Final := (1-(s0_CHA2DS2_VASc_Final)^exp(CHA2DS2_VASc_Final - (score_avgbeta_CHA2DS2_VASc_Final)))*100]; summary(mydata$pred_risk_CHA2DS2_VASc_Final)
mydata[,pred_risk_C2HEST_Final := (1-(s0_C2HEST_Final)^exp(C2HEST_Final - (score_avgbeta_C2HEST_Final)))*100]; summary(mydata$pred_risk_C2HEST_Final)

summary(mydata$EHR_AF_Final)
summary(mydata$CHARGE_AF_Final)
summary(mydata$CHA2DS2_VASc_Final)
summary(mydata$C2HEST_Final)

current_date = Sys.Date()

#######################################
### Calibration
#######################################

score_names = c("EHR_AF_Final", "CHARGE_AF_Final", "CHA2DS2_VASc_Final", "C2HEST_Final")
#score_names = c("CHARGE_AF_Final")

for (score_name in score_names)
{
  surmod <- with(mydata, Surv(num_days_to_AFib_or_to_cencor,is_AFib_combined))
  
  tmp_str_1 = as.formula(paste("surmod", score_name, sep = "~"))
  fit <- cph(tmp_str_1, data = mydata, surv = TRUE, time.inc = 5*365, u = 5*365, x = T, y = T)
  
  fit
  
  memory.limit(size = 456000)
  cal.hare <- calibrate(fit, u = 5*365, cmethod = 'hare', B = 1) # 200
  
  plot(cal.hare)
  plot(1 - cal.hare, xlab = "Predicted 5-year risk of AF", ylab = "Proportion with AF at 5 years", riskdist = TRUE, xlim = c(0, 0.8), ylim = c(0, 0.8), cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5)
  
  #Calibration estimate
  cal_charge_prs <- c(fit$coefficients[1], confint(fit, score_name)[1], confint(fit, score_name)[2])
  cal_charge_prs ## THIS IS THE CALIBRATION SLOPE, 95% LOWER BOUND, 95% UPPER BOUND
  
  #score_name = "CHARGE_AF_Final"
  tmp_file_name_1 = paste(score_name, current_date, sep = "_")
  file_name_1 = paste(tmp_file_name_1, "_Calibration.jpg", sep = "")
  
  dev.copy(jpeg, filename = file_name_1, width = 600, height = 600);
  dev.off()
  
  #######################################
  ###Risk plot
  #######################################
  
  tmp_str_2 = paste('pred_risk_', score_name, sep ='')
  quantile(mydata[[tmp_str_2]], probs = c(0, 0.33, 0.67, 1))
  
  mydata$Group <- ifelse(mydata[[tmp_str_2]] > 5, ">=5%", ifelse(mydata[[tmp_str_2]] < 2.5, "<2.5%", "2.5%-5%"))
  
  count(mydata$Group)
  
  mydata$Group <- factor(mydata$Group, levels = c(">=5%", "2.5%-5%", "<2.5%"))
  
  ci.af_v <- prodlim(Hist(num_days_to_AFib_or_to_cencor, is_AFib_combined) ~  Group, data = mydata)
  plot(ci.af_v, "cuminc", ylim = c(0, 0.15), xlab = "Follow-up (days)", ylab = "Cumulative AF risk (%)", col = c("red", "orange", "yellow"), cex.axis = 5, cex.main = 5, cex.lab = 5)
  
  tmp_file_name_2 = paste(score_name, current_date, sep = "_")
  file_name_2 = paste(tmp_file_name_2, "_Risk_Groups.jpg", sep = "")
  
  dev.copy(jpeg, filename = file_name_2, width = 1024, height = 768);
  dev.off()
  
}

#######################################
# Generating a density plot
#######################################

for (score_name in score_names)
{
  score_name = paste("pred_risk_", score_name, sep = "")
  
  score_name_trimmed = str_replace(score_name, '_Final', '');
  density_str_1 = "Predicted 5-year AF risk using "; density_str_2 = paste(density_str_1, score_name_trimmed, sep = ""); density_str_3 = paste(density_str_2, " score (%)", sep = "")
  # Create separate data frames for cases/controls
  incident_af <- mydata[is_AFib_combined == 1,]
  no_incident_af <- mydata[is_AFib_combined == 0,]
  x <- list(v1=incident_af[[score_name]],v2=no_incident_af[[score_name]])
  data <- melt(x)
  # Density of predicted risk distribution
  tmp_density_file_name_1 = paste(score_name, current_date, sep = "_Density_Plot_")
  tmp_density_file_name_2 = paste(tmp_density_file_name_1, '.png', sep = "")
  png(file = tmp_density_file_name_2,height=540,width=760)
  ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) +
    scale_x_continuous(breaks=seq(0,20,1),expand=c(0,0.1),limits=c(0,20)) +
    scale_y_continuous(breaks=seq(0,0.35,0.05),expand=c(0,0),limits=c(0,0.35)) +
    scale_fill_manual(values=c("#2b8cbe","#f03b20"), name='',labels=c('Incident AF', 'No incident AF')) +
    theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
          axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
          axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
    xlab(density_str_3) + ylab("Density")
  dev.off()
  ggsave(tmp_density_file_name_2)
}

## Histogram

mydata_only_patients_with_AFib_outcome = mydata[which(mydata$is_AFib_combined == 1),]
tmp1_for_hist = as.numeric(mydata_only_patients_with_AFib_outcome$num_days_to_AFib_or_to_cencor)

hist(tmp1_for_hist, breaks = 100, main = "AFib population", xlab = "Days from index date", ylab = "# AFib patients", xaxt = 'n')
axis(side = 1, at = c(0, 365, 730, 1095, 1460, 1825))

#######################################