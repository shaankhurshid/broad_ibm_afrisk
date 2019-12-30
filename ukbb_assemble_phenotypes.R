# Script for building initial UKBB dataset for clinical scores/GRS

# Dependencies
library(data.table)

# Load initial dataset
load("/Volumes/medpop_afib/skhurshid/PRS/ukbb_output.RData")

############################ MERGE INITIAL DATASETS ############################
## Create list of names
names <- list(af=af,cad=cad,ckd=ckd,dm=dm,hf=hf,hld=hld,htn=htn,hypothyroid=hypothyroid,is=is,
              mi=mi,pad=pad,stroke=stroke,vd=vd)

## Looping cbind
phenos <- data.frame(ID=af$ID)

for (i in 1:length(names)){
  phenos <- cbind(phenos,names[[i]][,2:7])
  colnames(phenos)[(length(phenos)-2):length(phenos)] <-
    paste(names(names)[i],colnames(phenos)[(length(phenos)-2):length(phenos)],sep='')
}

############################ SEX ############################
# Script for extracting baseline SBP from raw data

# Read in sex values
sex <- fread(file="ukb9222_f.31.tab")

# Convert to factor
sex[,sex:=ifelse(f.31.0.0==0,'female','male'),by=f.eid]

# Merge with parent
phenos$sex <- sex$sex

############################ TOBACCO PHENOTYPE ############################

# Read in raw ICD9 and ICD10 codes
icd9_main <- read.delim(file="ukb9222_f.41203.tab",sep="\t")
icd9_secondary <- read.delim(file="ukb9222_f.41205.tab",sep="\t")
icd10_main <- read.delim(file="ukb9222_f.41202.tab",sep="\t")
icd10_secondary <- read.delim(file="ukb9222_f.41204.tab",sep="\t")

# Define Codes
tobacco_codes_ICD9 <- as.factor(3051)
tobacco_codes_ICD10 <- as.factor(c("F170","F171","F172","F173","T652"))

# ICD 9
tobacco_index_ICD9_main <- apply(icd9_main,2,function(x) x %in% tobacco_codes_ICD9)
tobacco_index_ICD9_secondary <- apply(icd9_secondary,2,function(x) x %in% tobacco_codes_ICD9)

# ICD 10
tobacco_index_ICD10_main <- apply(icd10_main,2,function(x) x %in% tobacco_codes_ICD10)
tobacco_index_ICD10_secondary <- apply(icd10_secondary,2,function(x) x %in% tobacco_codes_ICD10)

# Get sums
tobacco_sum <- rowSums(tobacco_index_ICD9_main) + rowSums(tobacco_index_ICD9_secondary) + rowSums(tobacco_index_ICD10_main) + rowSums(tobacco_index_ICD10_secondary)

# Obtain ordinal variables
tobacco_ordinal <- read.delim(file='ukb9222_f.20116.tab')

# Consolidate sums with ordinals (2 = Current, 1 = Previous, 0 = Never)
tobacco <- cbind(tobacco_ordinal,tobacco_sum)
tobacco$tobacco_phenotype <- with(tobacco,
                                  ifelse(c((f.20116.0.0 == 2) | (tobacco_sum >= 1)),2,
                                         ifelse((f.20116.0.0 == 1),1,
                                                ifelse((f.20116.0.0 == 0),0,NA))))

# Reformat tobacco for merge
colnames(tobacco)[1] <- "ID"
tobacco <- tobacco[,c(1,6)]

# Merge with parent file
phenos <- cbind(phenos,tobacco["tobacco_phenotype"])

# Rename and factor
phenos$tobacco_phenotype[phenos$tobacco_phenotype==0] <- 'Never'
phenos$tobacco_phenotype[phenos$tobacco_phenotype==1] <- 'Previous'
phenos$tobacco_phenotype[phenos$tobacco_phenotype==2] <- 'Current'
phenos$tobacco_phenotype <- factor(phenos$tobacco_phenotype,levels=c('Current','Previous','Never'))

# Binary tobacco
phenos[,tobacco_binary := ifelse(tobacco_phenotype=="Current","Active","Non-Active")]

############################ RACE PHENOTYPE ############################
# Script for extracting race codes from raw data

# Read in ICD9 and ICD10 codes
race_raw <- read.delim(file="ukb9222_f.21000.tab",sep="\t")

# Make race phenotype
race <- race_raw[,c(1,2)]
race$race_phenotype <- with(race,
                            ifelse(c((!is.na(f.21000.0.0)) & c((f.21000.0.0 == 1) | (f.21000.0.0 == 1001) | (f.21000.0.0 == 1002)
                                                               | (f.21000.0.0 == 1003))),"White",
                                   ifelse(c((!is.na(f.21000.0.0)) & (f.21000.0.0 == 3001)),"Indian",
                                          ifelse(c((!is.na(f.21000.0.0)) & (f.21000.0.0 == 4001)),"Caribbean",
                                                 ifelse(c((!is.na(f.21000.0.0)) & (f.21000.0.0 == 4002)),"African",
                                                        ifelse(is.na(f.21000.0.0),NA,"Nonwhite"))))))

# Reformat race for merge
colnames(race)[1] <- "ID"
race <- race[,c(1,3)]

# Merge with parent file
phenos <- cbind(phenos,race["race_phenotype"])

# Rename and factor
phenos$race_phenotype <- factor(phenos$race_phenotype,levels=c('White','Nonwhite','Indian','Caribbean','African'))

# Race binary
phenos$race_binary <- with(phenos,
                           ifelse(c(!is.na(race_phenotype) & (race_phenotype=='White')),'White',
                                  ifelse(is.na(race_phenotype),NA,'Nonwhite')))

############################ SBP PHENOTYPE
### AUTO
# Script for extracting SBP from raw data

# Read in SBP values
sbp_raw <- fread(file="ukb9222_f.4080.tab")
# Baseline only
sbp_baseline_auto <- sbp_raw[,.(ID=f.eid, aSBP1 = f.4080.0.0, aSBP2 = f.4080.0.1)]

### MANUAL
# Read in values
mbp <- fread("/Volumes/medpop_afib/jhalford/UKBB_prs_etoh/data/manualbp.txt")
names(mbp) <- c("ID", "mSBP1", "mSBP2", "mDBP1", "mDBP2")
sbp_baseline_man <- mbp[,c('ID','mSBP1','mSBP2')]

### COMBINE
# Set keys for merge
setkey(sbp_baseline_auto,ID)
setkey(sbp_baseline_man,ID)

# Perform left join
sbp_baseline <- merge(sbp_baseline_auto,sbp_baseline_man,all.x=TRUE)
sbp_average <- sbp_baseline[,.(SBP1 = mean(c(aSBP1,aSBP2,mSBP1,mSBP2),na.rm=T)),by=ID]

# Merge with main dataset

# Set keys for merge
setkey(sbp_average,ID)
setkey(phenos,ID)

# Perform left join
phenos <- phenos[sbp_average]

############################ DBP PHENOTYPE
# Script for extracting DBP from raw data

# Read in DBP values
dbp_raw <- fread(file="ukb9222_f.4079.tab")

# Baseline only
dbp_baseline_auto <- dbp_raw[,.(ID=f.eid, aDBP1 = f.4079.0.0, aDBP2 = f.4079.0.1)]
dbp_baseline_man <- mbp[,c('ID','mDBP1','mDBP2')]

### COMBINE
# Set keys for merge
setkey(dbp_baseline_auto,ID)
setkey(dbp_baseline_man,ID)

# Perform left join
dbp_baseline <- merge(dbp_baseline_auto,dbp_baseline_man,all.x=TRUE)
dbp_average <- dbp_baseline[,.(DBP1 = mean(c(aDBP1,aDBP2,mDBP1,mDBP2),na.rm=T)),by=ID]

# Merge with main dataset

# Set keys for merge
setkey(dbp_average,ID)
setkey(phenos,ID)

# Perform left join
phenos <- phenos[dbp_average]

############################ Height PHENOTYPE
# Script for extracting height from raw data

# Read in height values
ht_raw <- fread(file="ukb9222_f.50.tab")

# Any height measurement available (if multiple, take average)
ht <- ht_raw[,apply(.SD,1,mean,na.rm=T),by=f.eid]

# Merge with main dataset
phenos$ht <- ht$V1

############################ Weight PHENOTYPE
# Script for extracting weight from raw data

### Standard weight
# Read in weight values
wt_raw <- fread(file="ukb9222_f.21002.tab")

# Any weight measurement available (if multiple, take average)
wt <- wt_raw[,apply(.SD,1,mean,na.rm=T),by=f.eid]

# Merge with main dataset
phenos$wt_std <- wt$V1

### Impedance weight
# Read in weight values
wt_raw <- fread(file="ukb9222_f.23098.tab")

# Any weight measurement available (if multiple, take average)
wt <- wt_raw[,apply(.SD,1,mean,na.rm=T),by=f.eid]

# Merge with main dataset
phenos$wt_imp <- wt$V1

# Weight combined (if multiple, take average)
phenos <- data.table(phenos)
wt_all <- phenos[,apply(.SD,1,mean,na.rm=T),by=ID,.SDcols=c("wt_std","wt_imp")]

# Merge with main dataset
phenos$wt_all <- wt_all$V1

############################ BP MEDS PHENOTYPE
## Medication field 1: 6153
# Read in medication values
bpmeds_raw <- fread(file="ukb9222_f.6153.tab")

# Mark for BP medication at baseline
## Dummy function that searches for 2s (2 = Blood Pressure Medication)
find_bp <- function(x){
  value <- ifelse(c(c(!is.na(x[1]) & (x[1] == 2)) | c(!is.na(x[2]) & (x[2] == 2))
                    | c(!is.na(x[3]) & (x[3] == 2)) | c(!is.na(x[4]) & (x[4] == 2))),1,0)
}

bpmeds <- bpmeds_raw[,apply(.SD,1,find_bp),by=f.eid]
phenos$bp_med1 <- bpmeds$V1

## Medication field 2: 6177
# Read in medication values
bpmeds_raw <- fread(file="ukb9222_f.6177.tab")

# Mark for BP medication at baseline
## Dummy function that searches for 2s (2 = Blood Pressure Medication)
find_bp <- function(x){
  value <- ifelse(c(c(!is.na(x[1]) & (x[1] == 2)) | c(!is.na(x[2]) & (x[2] == 2))
                    | c(!is.na(x[3]) & (x[3] == 2))),1,0)
}

bpmeds2 <- bpmeds_raw[,apply(.SD,1,find_bp),by=f.eid]
phenos$bp_med2 <- bpmeds2$V1

## Combined BP medication variable
phenos$bp_med_combined <- with(phenos,
                               ifelse(c((bp_med1==1) | (bp_med2==1)),1,0))

############################ AF AND STROKE TIME
# AF FU time
phenos[,af_fu_time := afagecensor - afagevisit0, by=ID]

# Stroke FU time
phenos[,stroke_fu_time := strokeagecensor - strokeagevisit0, by=ID]

############################ Remove NAs from prevalent phenotypes
names <- c('Prev_Heart_Failure_V2','Prev_Hypertension','Prev_Diabetes_All','Prev_Stroke','Prev_Coronary_Artery_Disease',
           'Prev_Peripheral_vascular_disease','Prev_Myocardial_Infarction','Prev_Hypercholesterolemia',
           'Prev_Valve_Disease','Prev_Chronic_kidney_disease','Prev_Hypothyroidism')
for (j in names){set(phenos,i=which(is.na(phenos[[j]])),j=j,value=0)}

############################ CHADS Score
## At visit zero
phenos[,':='(chads_visit0 = (Prev_Heart_Failure_V2==1) + (Prev_Hypertension==1)
                  + ifelse(afagevisit0 >= 75,2,ifelse(afagevisit0 >= 65,1,0))
                  + (Prev_Diabetes_All==1) + (sex=='female') + (Prev_Stroke==1)*2
                  + ifelse(c((Prev_Coronary_Artery_Disease==1) | (Prev_Peripheral_vascular_disease==1) | (Prev_Myocardial_Infarction==1)),1,0))]

## At time of stroke
# Create prevalent at time of stroke variables
phenos[Incd_Stroke==1,':='(age_at_stroke = strokeagecensor,
        prev_hf_atstroke = ifelse(is.na(Incd_Heart_Failure_V2),1,
                                  ifelse(c(Incd_Heart_Failure_V2==1 & (hfagecensor < strokeagecensor)),1,0)),
        prev_htn_atstroke = ifelse(is.na(Incd_Hypertension),1,
                                   ifelse(c(Incd_Hypertension==1 & (htnagecensor < strokeagecensor)),1,0)),
        prev_dm_atstroke = ifelse(is.na(Incd_Diabetes_All),1,
                                  ifelse(c(Incd_Diabetes_All==1 & (dmagecensor < strokeagecensor)),1,0)),
        prev_cad_atstroke = ifelse(is.na(Incd_Coronary_Artery_Disease),1,
                                   ifelse(c(Incd_Coronary_Artery_Disease==1 & (cadagecensor < strokeagecensor)),1,0)),
        prev_mi_atstroke = ifelse(is.na(Incd_Myocardial_Infarction),1,
                                  ifelse(c(Incd_Myocardial_Infarction==1 & (miagecensor < strokeagecensor)),1,0)),
        prev_pad_atstroke = ifelse(is.na(Incd_Peripheral_vascular_disease),1,
                                   ifelse(c(Incd_Peripheral_vascular_disease==1 & (padagecensor < strokeagecensor)),1,0)))]

# Apply variables to generate score
phenos[Incd_Stroke==1,':='(chads_at_stroke = (prev_hf_atstroke==1) + (prev_htn_atstroke==1)
                  + ifelse(age_at_stroke >= 75,2,ifelse(age_at_stroke >= 65,1,0))
                  + (prev_dm_atstroke==1) + (sex=='female') +
                  + ifelse(c((prev_cad_atstroke==1) | (prev_pad_atstroke==1) | (prev_mi_atstroke==1)),1,0))]

## At time of AF
# Create prevalent at time of AF variables
phenos[Incd_Atrial_fibrillation_or_flutter==1,':='(age_at_AF = afagecensor,
                           prev_hf_atAF = ifelse(is.na(Incd_Heart_Failure_V2),1,
                                                     ifelse(c(Incd_Heart_Failure_V2==1 & (hfagecensor < afagecensor)),1,0)),
                           prev_htn_atAF = ifelse(is.na(Incd_Hypertension),1,
                                                      ifelse(c(Incd_Hypertension==1 & (htnagecensor < afagecensor)),1,0)),
                           prev_dm_atAF = ifelse(is.na(Incd_Diabetes_All),1,
                                                     ifelse(c(Incd_Diabetes_All==1 & (dmagecensor < afagecensor)),1,0)),
                           prev_cad_atAF = ifelse(is.na(Incd_Coronary_Artery_Disease),1,
                                                      ifelse(c(Incd_Coronary_Artery_Disease==1 & (cadagecensor < afagecensor)),1,0)),
                           prev_mi_atAF = ifelse(is.na(Incd_Myocardial_Infarction),1,
                                                     ifelse(c(Incd_Myocardial_Infarction==1 & (miagecensor < afagecensor)),1,0)),
                           prev_pad_atAF = ifelse(is.na(Incd_Peripheral_vascular_disease),1,
                                                      ifelse(c(Incd_Peripheral_vascular_disease==1 & (padagecensor < afagecensor)),1,0)))]

# Apply variables to generate score
phenos[Incd_Atrial_fibrillation_or_flutter==1,':='(chads_at_AF = (prev_hf_atAF==1) + (prev_htn_atAF==1)
                           + ifelse(age_at_stroke >= 75,2,ifelse(age_at_stroke >= 65,1,0))
                           + (prev_dm_atAF==1) + (sex=='female') +
                             + ifelse(c((prev_cad_atAF==1) | (prev_pad_atAF==1) | (prev_mi_atAF==1)),1,0))]

############################ CHADS Score
# Binary "merits anticoagulation" cutoff variable
## Baseline
phenos[,':='(needs_oac_visit0 =
                                      ifelse(c((sex=='female') & (chads_visit0 >= 3)),1,
                                           ifelse(c((sex=='male') & (chads_visit0 >= 2)),1,0)))]

## At stroke
phenos[Incd_Stroke==1,':='(needs_oac_at_stroke =
                                               ifelse(c((sex=='female') & (chads_at_stroke >= 3)),1,
                                             ifelse(c((sex=='male') & (chads_at_stroke >= 2)),1,0)))]

## At AF
phenos[Incd_Atrial_fibrillation_or_flutter==1,':='(needs_oac_at_af =
                                    ifelse(c((sex=='female') & (chads_at_stroke >= 3)),1,
                                    ifelse(c((sex=='male') & (chads_at_stroke >= 2)),1,0)))]

############################ Age cutoff
## Age 65
phenos[,':='(age65_visit0 =
                    ifelse(afagevisit0 >= 65,1,0))]

############################ AF RISK GROUPS
phenos[,':='(charge_pred_af_above5 = ifelse(charge.pred5_cal >= 5,1,0),
               ehraf_pred_af_above5 = ifelse(ehraf.pred5_cal >= 5,1,0))]

phenos[,':='(charge_pred_af_above13 = ifelse(charge.pred5_cal >= 12.9,1,0),
                  ehraf_pred_af_above13 = ifelse(ehraf.pred5_cal >= 12.9,1,0))]

phenos[,':='(charge_pred_af_above3 = ifelse(charge.pred5_cal >= 3.26,1,0),
                  ehraf_pred_af_above3 = ifelse(ehraf.pred5_cal >= 3.26,1,0))]

phenos[,':='(charge_prs_pred_af_above5 = ifelse(charge_prs_pred5 >= 5,1,0),
                  ehraf_prs_pred_af_above5 = ifelse(ehraf_prs_pred5 >= 5,1,0))]

phenos[,':='(charge_prs_pred_af_above3 = ifelse(charge_prs_pred5 >= 3.26,1,0),
                  ehraf_prs_pred_af_above3 = ifelse(ehraf_prs_pred5 >= 3.26,1,0))]

############################ AF 5-YR PHENOTYPE
phenos[,incd_af_5y := as.double(ifelse(c((Incd_Atrial_fibrillation_or_flutter==1) & (af_fu_time<=5)),1,0)), by=ID]
phenos[,incd_af_5y.t := as.double(ifelse(af_fu_time > 5,5,af_fu_time)), by=ID]

phenos[,incd_af_5y := as.double(ifelse(c((Incd_Atrial_fibrillation_or_flutter==1) & (af_fu_time<=5)),1,0)), by=ID]
phenos[,incd_af_5y.t := as.double(ifelse(af_fu_time > 5,5,af_fu_time)), by=ID]

############################ Stroke 5-YR PHENOTYPE
### All stroke
pg_complete_noprevaf[,stroke_5y := as.double(ifelse(c((Incd_Stroke==1) & (stroke_fu_time<=5)),1,0)), by=ID]
pg_complete_noprevaf[,stroke_5y.t := as.double(ifelse(stroke_fu_time > 5,5,stroke_fu_time)), by=ID]

pg_complete_noprevafstroke[,stroke_5y := as.double(ifelse(c((Incd_Stroke==1) & (stroke_fu_time<=5)),1,0)), by=ID]
pg_complete_noprevafstroke[,stroke_5y.t := as.double(ifelse(stroke_fu_time > 5,5,stroke_fu_time)), by=ID]

### Ischemic stroke
pg_complete_noprevaf[,istroke_fu_time := isagecensor - isagevisit0]
pg_complete_noprevafstroke[,istroke_fu_time := isagecensor - isagevisit0]

pg_complete_noprevaf[,istroke_5y := as.double(ifelse(c((Incd_Ischemic_stroke==1) & (istroke_fu_time<=5)),1,0)), by=ID]
pg_complete_noprevaf[,istroke_5y.t := as.double(ifelse(istroke_fu_time > 5,5,istroke_fu_time)), by=ID]

pg_complete_noprevafstroke[,istroke_5y := as.double(ifelse(c((Incd_Ischemic_stroke==1) & (istroke_fu_time<=5)),1,0)), by=ID]
pg_complete_noprevafstroke[,istroke_5y.t := as.double(ifelse(istroke_fu_time > 5,5,istroke_fu_time)), by=ID]

############################ INCD STROKE BEFORE AF PHENOTYPE
### Before AF (any)
phenos[,stroke_before_af := as.double(ifelse(c((Incd_Stroke==1) & (afagecensor > strokeagecensor)),1,0))]
phenos[,stroke_before_af.t := ifelse(is.na(stroke_before_af),afagecensor - afagevisit0,
                                          ifelse(stroke_before_af==1,strokeagecensor - strokeagevisit0,
                                             ifelse(is.na(Incd_Atrial_fibrillation_or_flutter),stroke_fu_time,
                                                 ifelse(Incd_Atrial_fibrillation_or_flutter==1,afagecensor - afagevisit0,stroke_fu_time))))]

## Before AF (90 days)
pg_complete_noprevafstroke[,stroke_before_af90 := ifelse(c((stroke_before_af==1) & ((afagecensor - strokeagecensor) <= 90/365.25)),1,0)]
pg_complete_noprevafstroke[,stroke_before_af90.t := ifelse(stroke_before_af90==1,stroke_fu_time, # stroke before AF 90 - censor at outcome
                                     ifelse(c(Incd_Stroke==1 & Incd_Atrial_fibrillation_or_flutter==1),pmin((afagecensor - afagevisit0),(strokeagecensor - strokeagevisit0 + 90/365.25)), # stroke and AF but not outcome - censor at earliest of stroke + 90 days or AF
                                            ifelse(c(Incd_Stroke==1 & Incd_Atrial_fibrillation_or_flutter==0),(strokeagecensor - strokeagevisit0 + 90/365.25), # stroke but no AF - censor at stroke + 90 days
                                                   ifelse(c(Incd_Stroke==0 & Incd_Atrial_fibrillation_or_flutter==1),af_fu_time,af_fu_time))))] # AF but no stroke / No AF or stroke, censor as per AF

## Before AF (90 days, censored at 5 years)
pg_complete_noprevafstroke[,stroke_before_af90_5 := ifelse(c((stroke_before_af==1) & (stroke_5y==1) & ((afagecensor - strokeagecensor) <= 90/365.25)),1,0)]
pg_complete_noprevafstroke[,stroke_before_af90_5.t := ifelse(stroke_before_af90_5==1,stroke_5y.t, # stroke before AF 90 - censor at outcome
                                                           ifelse(c(Incd_Stroke==1 & Incd_Atrial_fibrillation_or_flutter==1),pmin((afagecensor - afagevisit0),(strokeagecensor - strokeagevisit0 + 90/365.25),5), # stroke and AF but not outcome - censor at earliest of stroke + 90 days or AF or 5 years
                                                                  ifelse(c(Incd_Stroke==1 & Incd_Atrial_fibrillation_or_flutter==0),pmin((strokeagecensor - strokeagevisit0 + 90/365.25),5), # stroke but no AF - censor at stroke + 90 days or 5 years
                                                                         ifelse(c(Incd_Stroke==0 & Incd_Atrial_fibrillation_or_flutter==1),pmin(af_fu_time,5),pmin(af_fu_time,5)))))] # AF but no stroke / No AF or stroke, censor as per AF

############################ No prevalent or incident AF dataset
phenos_noaf <- phenos[c((Prev_Atrial_fibrillation_or_flutter==0) & (Incd_Atrial_fibrillation_or_flutter==0)),]

############################ No prevalent AF dataset
phenos_noprevaf <- phenos[Prev_Atrial_fibrillation_or_flutter==0,]

############################ No prevalent or incident AF dataset
phenos_noaf <- phenos[c((Prev_Atrial_fibrillation_or_flutter==0) & (Incd_Atrial_fibrillation_or_flutter==0)),]

# Save output
save(phenos,file='phenos_080919.RData')
save(phenos_noprevaf,file='phenos_genos_noprevaf_080919.RData')
save(phenos_noaf,file='phenos_genos_noaf_080919.RData')
