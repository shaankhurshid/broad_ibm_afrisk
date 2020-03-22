rm(list=ls()); memory.limit(size = 456000); library(ROCR); library(ggplot2); library(pROC); library(stringr); library(plyr); library(prodlim); library('rms'); library(data.table); library(survival); library(reshape2); library(rms)

##############
### Parameters
#FILE = 'Random_Sample_AFib_Final_SUPERMART_122_Date_2020-02-25.csv'
FILE = 'Second_Data_Selection_Approach_AFib_Final_SUPERMART_122_Date_2020-02-17.csv'
DIRECTORY = "Y:/ibm-rwe_global/URI/Broad/Validation_Tool/output/"; setwd(DIRECTORY)
OUTCOME = 'is_AFib_combined'
NUM_DAYS_TO_OUTCOME = 'num_days_to_AFib_or_to_cencor'
SCORES = c('EHR_AF_Final', 'CHA2DS2_VASc_Final', 'CHARGE_AF_Final', 'C2HEST_Final')
#SCORES = c('EHR_AF_Final')
NUM_DAYS_FOLLOW_UP = 1825
B_calibrate_parameter = 10
AGE_MIN_YEARS = 45; AGE_MAX_YEARS = 90; AGE_STEP_YEARS = 5
##############

mydata <- read.csv(file = FILE, header = TRUE, sep = ','); setDF(mydata); #colnames(mydata)
AFib_incidence_in_entire_population = round(100 * nrow(mydata[which(mydata$is_AFib_combined == 1),]) / nrow(mydata), 2);

df_scores_final <- NULL

#Select population
for(is_population_restricted_by_age in seq(0, 1))
{
  if (is_population_restricted_by_age == 0)
    current_AGE_STEP_YEARS = AGE_MAX_YEARS - AGE_MIN_YEARS
  else
    current_AGE_STEP_YEARS = AGE_STEP_YEARS
  
  #print(current_AGE_STEP_YEARS)

for (gender_item in seq(0, 2))
{
  if (gender_item == 0)
    current_gender = 'Female-only'
  if (gender_item == 1)
    current_gender = 'Male-only'
  if (gender_item == 2)
    current_gender = 'All'
  
for (age_item in seq(AGE_MIN_YEARS, AGE_MAX_YEARS - current_AGE_STEP_YEARS, current_AGE_STEP_YEARS))
{
  if (age_item == AGE_MAX_YEARS - current_AGE_STEP_YEARS)
  {
    min_age_for_final_df = age_item; max_age_for_final_df = age_item + current_AGE_STEP_YEARS
    if (gender_item == 0 | gender_item == 1)
    {
    tmp_mydata = subset(mydata, age >= age_item & age <= age_item + current_AGE_STEP_YEARS & is_male == gender_item)
    Population_type_str = paste(paste(paste(paste(paste(paste('n = ', nrow(tmp_mydata), sep = ''), ' ', sep = ''), current_gender, sep = ''), ' Age ', age_item, sep = ''), ' to ', sep = ''), age_item + current_AGE_STEP_YEARS, sep ='')
    }
    else
    {
      tmp_mydata = subset(mydata, age >= age_item & age <= age_item + current_AGE_STEP_YEARS)
      Population_type_str = paste(paste(paste(paste(paste(paste('n = ', nrow(tmp_mydata), sep = ''), ' ', sep = ''), current_gender, sep = ''), ' Age ', age_item, sep = ''), ' to ', sep = ''), age_item + current_AGE_STEP_YEARS, sep ='')
    }
  }
  else
  {
    min_age_for_final_df = age_item; max_age_for_final_df = age_item + current_AGE_STEP_YEARS - 1
    if (gender_item == 0 | gender_item == 1)
    {
    tmp_mydata = subset(mydata, age >= age_item & age < age_item + current_AGE_STEP_YEARS & is_male == gender_item)
    Population_type_str = paste(paste(paste(paste(paste(paste('n=', nrow(tmp_mydata), sep = ''), ' ', sep = ''), current_gender, sep = ''), ' Age ', age_item, sep = ''), ' to ', sep = ''), age_item + current_AGE_STEP_YEARS - 1, sep ='')
    }
    else
    {
      tmp_mydata = subset(mydata, age >= age_item & age < age_item + current_AGE_STEP_YEARS)
      Population_type_str = paste(paste(paste(paste(paste(paste('n = ', nrow(tmp_mydata), sep = ''), ' ', sep = ''), current_gender, sep = ''), ' Age ', age_item, sep = ''), ' to ', sep = ''), age_item + current_AGE_STEP_YEARS - 1, sep ='')
    }
  }
  
  print(Population_type_str)
  
  for (SCORE_NAME in SCORES)
  {
    names(tmp_mydata)[names(tmp_mydata) == NUM_DAYS_TO_OUTCOME] <- 'time'; names(tmp_mydata)[names(tmp_mydata) == OUTCOME] <- 'status'; names(tmp_mydata)[names(tmp_mydata) == SCORE_NAME] <- 'score'
    
    max_follow_up_days_tmp_mydata = max(tmp_mydata$time)
    
    #Calculate prevalence
    AFib_incidence_in_subgroup = round(100 * nrow(tmp_mydata[which(tmp_mydata$status == 1),]) / nrow(tmp_mydata), 2);
    
    #Calculate calibration slope
    survmod <- with(tmp_mydata, Surv(time, status))
    fit<-cph(survmod~score, data = tmp_mydata, surv=TRUE, time.inc = max_follow_up_days_tmp_mydata, u = max_follow_up_days_tmp_mydata, x = T, y = T)
    fit2 <- coxph(survmod~score,data = tmp_mydata)
    cal_score <- c(fit2$coefficients[1], confint(fit2,"score")[1], confint(fit2,"score")[2]) # beta and 95%CI corresponds to calibration slope.
    calibration_results = round(cal_score, 3)
    
    #Calculate CI
    summary_cox_object = summary(fit2)
    concordance_index_low_conf = round(summary_cox_object$concordance[1] - 1.96 * summary_cox_object$concordance[2], 3)
    concordance_index_value = round(summary_cox_object$concordance[1], 3)
    concordance_index_high_conf = round(summary_cox_object$concordance[1] + 1.96 * summary_cox_object$concordance[2], 3)
    
    ##Add values to summary data frame
    age_range_str = paste(min_age_for_final_df, max_age_for_final_df, sep = '-')
    df_scores_final <-rbind(df_scores_final, data.frame(n = nrow(tmp_mydata), Population_Type = current_gender, Min_Age = min_age_for_final_df, Max_Age = max_age_for_final_df, Age_Range = age_range_str, Score_Name = SCORE_NAME, AFib_Incidence_in_Subgroup = AFib_incidence_in_subgroup, AFib_Incidence_in_Entire_Population = AFib_incidence_in_entire_population, Calibration_Slope_Low_Conf = calibration_results[2], Calibration_Slope_Value = calibration_results[1], Calibration_Slope_High_Conf = calibration_results[3], Concordance_Index_Low_Conf = concordance_index_low_conf, Concordance_Index_Value = concordance_index_value, Concordance_Index_High_Conf = concordance_index_high_conf))
    
    # plot calibration
    cal.score <- 999
    while (cal.score == 999)
    {
    cal.score = tryCatch(
      {calibrate(fit, u = max_follow_up_days_tmp_mydata, cmethod='hare', B = B_calibrate_parameter)},
      error = function(a)
      {
        print(999)
      }
    )
    }
    
    SCORE_NAME = str_replace(SCORE_NAME, '_Final', '')
    file_name = paste(paste(Population_type_str, SCORE_NAME, sep = "_"), '.pdf', sep = '')
    
    pdf(file_name,height=3,width=3,pointsize=3); par(oma=c(1,1,1,1)); par(oma=c(1,1,1,1));col2=paste(rgb(252,146,114,maxColorValue=255),sep="")
    plot(cal.score,scat1d.opts=list(frac=0.1,side=1),xlim=c(1,0.2),ylim=c(1,0.2),xaxt="n",yaxt="n", xlab="Predicted 5-year risk of AF", ylab="Proportion with AF at 5 years", subtitles=F, par.corrected=list(col=col2, lty=1, lwd=2), bty='n', cex.lab=1.25)
    axis(1,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0)); axis(2,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0), las=1)
    mtext(text=SCORE_NAME, side=3, line=-0.5, adj=0.5, cex=1.2);
    mtext(text=Population_type_str, side=3, line=-1.6, adj=0.5, cex=1.2);
    mtext(text="Calibration slope", side=3, line=-2.5, at=0.85, cex=0.9); cal_rounded = round(cal_score, 3); calibration_slope_text = paste(paste(paste(paste(cal_rounded[1], " (", sep = ""), cal_rounded[2], "-", sep = ""), cal_rounded[3], sep = ""), ")", sep = "")
    mtext(text = calibration_slope_text, side=3, line=-3.5, at=0.85, cex=0.9);
    
    mtext(text="Concordance index", side=3, line=-4.8, at=0.85, cex=0.9); concordance_index_text = paste(paste(paste(paste(concordance_index_low_conf, " (", sep = ""), concordance_index_value, "-", sep = ""), concordance_index_high_conf, sep = ""), ")", sep = "")
    mtext(text = concordance_index_text, side=3, line=-5.8, at=0.85, cex=0.9);
    
    legend(0.5,0.8,c('Optimal','Observed','Optimism-corrected'),col=c('darkgray','black',col2),lty=1,lwd=1.5,bty='n',cex=0.75)
    segments(0,0,0.9,0.9,col='darkgray',lty=1);dev.off()
    
    tmp_mydata$score <- NULL
  }
  ###
  }
}
}

write.table(df_scores_final, file = "Calibration_results.csv", quote = FALSE, row.names = FALSE, sep = '\t')

