# Script for generating clinical scores for UKBB dataset
################################################# CREATE SCORE VARIABLES
pg[,':='(age_5 = afagevisit0/5,age_10 = afagevisit0/10, age_10_sq = (afagevisit0/10)**2,
                   ht_10 = ht/10, wt_15 = wt_all/15, sbp_20 = sbp0/20, dbp_10 = dbp0/10,
                   ht_10_sq = (ht/10)**2, wt_15_sq = (wt_all/15)**2,
                   dbp_less80 = as.double(dbp0<80)),by=ID]

################################################# CHARGE-AF
# Complete case
charge_variables <- as.character(c('afagevisit0','ht','wt_all','sbp0','dbp0',
                                   'race_binary','tobacco_phenotype',
                                   'bp_med_combined','Prev_Diabetes_All','Prev_Heart_Failure_V2',
                                   'Prev_Myocardial_Infarction'))

pg[,charge_complete := ifelse(c(is.na(afagevisit0) | is.na(ht) | is.na (wt_all) |
                                        is.na(sbp0) | is.na(dbp0) |
                                        is.na(race_binary) | is.na(tobacco_phenotype)
                                      | is.na(bp_med_combined) | is.na(Prev_Diabetes_All) |
                                        is.na(Prev_Heart_Failure_V2) | is.na(Prev_Myocardial_Infarction)),0,1),by=ID]

# CHARGE-AF
pg[charge_complete==1,charge := afagevisit0/5*0.508 + white*0.465 + ht/10*0.248 + wt_all/15*0.115 + sbp0/20*0.197
      +dbp0/10*(-0.101)+(tobacco_phenotype=="Current")*0.359+bp_med_combined*0.349+Prev_Diabetes_All*0.237
      +Prev_Heart_Failure_V2*0.701+Prev_Myocardial_Infarction*0.496,by=ID]

# Standardized CHARGE-AF
charge_mean <- pg[charge_complete==1,mean(charge)]
charge_sd <- pg[charge_complete==1,sd(charge)]
pg[charge_complete==1,charge_std := (charge-(charge_mean))/(charge_sd),by=ID]

# Predicted Risk for CHARGE-AF
pg[charge_complete==1,charge.pred5 := (1-0.9718412736^exp(charge-12.5815600))*100,by=ID]

################################################# EHR-AF
# Complete case
pg[,ehr_complete := ifelse(c(is.na(sex) | is.na(afagevisit0) | is.na(ht) | is.na (wt_all) |
                                      is.na(dbp0) | is.na(Prev_Hypertension) | is.na(race_binary) | is.na(tobacco_phenotype) |
                                      is.na(Prev_Hypercholesterolemia) | is.na(bp_med_combined) | is.na(Prev_Diabetes_All) |
                                      is.na(Prev_Heart_Failure_V2) | is.na(Prev_Coronary_Artery_Disease) |
                                      is.na(Prev_Valve_Disease) | is.na(Prev_Stroke) | is.na(Prev_Peripheral_vascular_disease) |
                                      is.na(Prev_Chronic_kidney_disease) | is.na(Prev_Hypothyroidism)),0,1),by=ID]

# EHR-AF
pg[,ehraf := (sex=="female")*(-0.137) + afagevisit0/10*1.494 + (afagevisit0/10)**2*(-0.048)
       + (race_binary=="White")*(-0.208) + (tobacco_phenotype=="Current")*0.152 + ht/10*(-0.231) + (ht/10)**2*(0.012)
       + wt_all/15*(-0.050) + (wt_all/15)**2*(0.021) + (dbp0<=80)*(-0.104) + Prev_Hypertension*0.106
       + Prev_Hypercholesterolemia*(-0.156) + Prev_Heart_Failure_V2*0.563 + Prev_Coronary_Artery_Disease*0.210
       + Prev_Valve_Disease*0.487 + Prev_Stroke*0.132 + Prev_Peripheral_vascular_disease*0.126 + Prev_Chronic_kidney_disease*0.279
       + Prev_Hypothyroidism*(-0.138)]

# Predicted Risk for EHR-AF
pg[ehr_complete==1,ehraf.pred5 := (1-0.9712209^exp(ehraf-6.728))*100,by=ID]

# Standardized EHR-AF
ehr_mean <- pg[ehr_complete==1,mean(ehraf)]
ehr_sd <- pg[ehr_complete==1,sd(ehraf)]
pg[ehr_complete==1,ehraf_std := (ehraf-(ehr_mean))/(ehr_sd),by=ID]

# Tertiles
## CHARGE-AF
pg_complete[,':='(charge_tertile = ifelse(charge < quantile(charge,probs=c(0.333,0.667))[1],
                                          "Low",ifelse(charge < quantile(charge,probs=c(0.333,0.667))[2],
                                                       "Medium","High")),
                  ehraf_tertile = ifelse(ehraf < quantile(ehraf,probs=c(0.333,0.667))[1],
                                         "Low",ifelse(ehraf < quantile(ehraf,probs=c(0.333,0.667))[2],
                                                      "Medium","High")),
                  prs_tertile = ifelse(prs < quantile(prs,probs=c(0.333,0.667))[1],
                                         "Low",ifelse(prs < quantile(prs,probs=c(0.333,0.667))[2],
                                                      "Medium","High")))]

################################################# 'STANDARDIZED PRS'
pg_complete_noprevaf[,prs_norm := prs/1126]

################################################# SIMPLE PRS
pg_complete_noprevaf[,simple_prs := charge + prs]

################################################# AGE
# Age score
pg_complete[,age_score := afagevisit0*0.112]

# Standardized age score
age_mean <- pg_complete[,mean(afagevisit0)]
age_sd <- pg_complete[,sd(afagevisit0)]
pg_complete[,age_std := (afagevisit0-(age_mean))/(age_sd),by=ID]

# Predicted Risk for CHARGE-AF
pg[charge_complete==1,charge.pred5 := (1-0.9718412736^exp(charge-12.5815600))*100,by=ID]
