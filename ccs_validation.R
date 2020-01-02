# Script to validate Explorys findings in PHS legacy cohort

# Dependencies
library(data.table)
library(stringr)
library(survival)
library(plyr)

###################################################### STEP 1: CLEAN LEGACY DATA
# Load legacy derivation
load(file='/data/arrhythmia/skhurshid/ehr_af/ds_021219.RData'); setDT(ds)

# Load legacy validation
load(file='/data/arrhythmia/skhurshid/ehr_af/vs_021219.RData'); setDT(vs)

# Reduce each only to necessary variables to reduce overhead
ds <- ds[,c('EMPI','start_fu','af_5y_sal','af_5y_sal.t','score','pred5')]
vs <- vs[,c('EMPI','start_fu','af_5y_sal','af_5y_sal.t','score','pred5')]

# Simple bind
all <- rbind(ds,vs)

###################################################### STEP 2: PROCESS CCS DATA
# Load and bind files
data <- list()
for (i in 1:length(list.files('/data/arrhythmia/skhurshid/ccs2'))){
  load(paste0('/data/arrhythmia/skhurshid/ccs2/',list.files('/data/arrhythmia/skhurshid/ccs2')[i]))
  data[[i]] <- s0[c(str_detect(names(s0),'._fin$') | str_detect(names(s0),'._fin_d$') | grepl('EMPI',names(s0)))]
  data[[i]] <- data[[i]][data[[i]][,'EMPI'] %in% all$EMPI,]
}

# Collapse
data2 <- do.call('rbind',data)

# Unique
ccs <- data2[!duplicated(data2$EMPI),]

# Join on legacy data
setDT(ccs)
setkey(ccs,'EMPI'); setkey(all,'EMPI')
ccs_validation <- ccs[all,nomatch=0]

# Throw away EMPI
ccs_validation <- ccs_validation[,!'EMPI']

# Save CCS data
save(ccs_validation,file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/ccs_validation.RData')

############################################################################################################
# STEPS 3-5 ARE CCS HIERARCHY SPECIFIC. YOU WILL NEED TO REPEAT FOR EACH HIERARCHICAL LEVEL #
# Just replace ccs1 with ccs2 or ccs3, etc
############################################################################################################

###################################################### STEP 3: Create prevalent variables based on start_fu
# Gather CCS1 related variables and park in list of names
names_ccs1 <- names(ccs_validation)[c(str_detect(names(ccs_validation),'ccs1') & 
                                        str_detect(names(ccs_validation),'_fin$'))]

str_sub(names_ccs1,-4,-1) <- '_baseline'

# Loop over the list to create baseline variables
## Set DF
setDF(ccs_validation)

## Create list of names for input
input_names <- names_ccs1
str_sub(input_names,-9,-1) <- ''

## Loop over the list (baseline = 1 if CCS condition occurred prior to start_fu)
for (i in 1:length(names_ccs1)){
  ccs_validation[,names_ccs1[i]] <- ifelse(c((ccs_validation[,paste0(input_names[i],'_fin')] == 1) &
                                             (ccs_validation[,paste0(input_names[i],'_fin_d')] <= ccs_validation[,'start_fu'])),1,0)
}

## Set back to DT
setDT(ccs_validation)

###################################################### STEP 4: Concordance analyses
# Summary stats for CCS1
## Condition counts
disease_counts <- ccs_validation[,lapply(.SD,sum),.SDcols=names_ccs1]

## AF counts
setDF(ccs_validation)
af_cts <- rep(NA,length(names_ccs1))

for (i in 1:length(names_ccs1)){
  af_cts[i] <- sum(ccs_validation[ccs_validation[,names_ccs1[i]]==1,]$af_5y_sal)
}

## Concordance
c <- rep(NA,length(names_ccs1))

for (i in 1:length(names_ccs1)){
  c[i] <- summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=ccs_validation[ccs_validation[,names_ccs1[i]]==1,]))$concordance[1]
}

## Concordance SE (requires bootstrap)
# Bootstrap function
boot <- function(time,status,response,data,runs,size=nrow(data)){
  out <- rep(NA,times=runs)
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE),]
    out[i] <- summary(coxph(Surv(sample[,time],sample[,status]) ~ sample[,response],data=sample))$concordance[1]
  }
  return(out)
}

## Condense output
output <- data.frame(ccs1=names_ccs1,n_disease=unlist(disease_counts),incd_af=af_cts,concordance=c,
                     concordance_se=se)

output$concordance_lb <- output$concordance-1.96*output$concordance_se
output$concordance_ub <- output$concordance+1.96*output$concordance_se

## Save out
write.csv(output,file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/ccs1_output.csv',row.names=F)

###################################################### STEP 5: Plotting
# Load output file with hand-curated disease names (OPTIONAL)
output <- read.csv('/data/arrhythmia/skhurshid/broad_ibm_afrisk/ccs1_output.csv')

# Sort data by concordance (low to high)
output <- output[order(output$concordance),]

png('/data/arrhythmia/skhurshid/broad_ibm_afrisk/concordance_ccs1.png',height=500,width=800)
par(oma=c(1,1,1,1))
par(mar=c(3,32,1,1))
plot(x=output$concordance,y=seq(1,18,1),xlim=c(0.700,1.00),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col='#6baed6',cex=1)
axis(1,at=seq(0.700,1.00,0.025),cex=1.4)
axis(2,at=1:18,labels=FALSE,cex=1.4)
mtext('CCS class (level 1)',side=2,line=31.5,cex=1.3)
mtext('Concordance',side=1,line=2.5,cex=1.3)
segments(output$concordance_lb,1:18,output$concordance_ub,1:18,col='#6baed6',lwd=2.2)

plot_names <- output$ccs1
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i],' (N=',output$n_disease[i],'; AF=',output$incd_af[i],')'),
        side=2,line=1,las=2,at=i)
}

text(output$concordance_ub+0.03,1:17,as.character(round(output$concordance,3)))
text(output$concordance_lb[18]-0.03,18,as.character(round(output$concordance[18],3)))
dev.off()

############################################################################################################
# CCS LEVEL 2 NOW
############################################################################################################

###################################################### STEP 3: Create prevalent variables based on start_fu
# Gather CCS1 related variables and park in list of names
names_ccs2 <- names(ccs_validation)[c(str_detect(names(ccs_validation),'ccs2') & 
                                        str_detect(names(ccs_validation),'_fin$'))]

str_sub(names_ccs2,-4,-1) <- '_baseline'

# Loop over the list to create baseline variables
## Set DF
setDF(ccs_validation)

## Create list of names for input
input_names <- names_ccs2
str_sub(input_names,-9,-1) <- ''

## Loop over the list (baseline = 1 if CCS condition occurred prior to start_fu)
for (i in 1:length(names_ccs2)){
  ccs_validation[,names_ccs2[i]] <- ifelse(c((ccs_validation[,paste0(input_names[i],'_fin')] == 1) &
                                               (ccs_validation[,paste0(input_names[i],'_fin_d')] <= ccs_validation[,'start_fu'])),1,0)
}

## Set back to DT
setDT(ccs_validation)

###################################################### STEP 4: Concordance analyses
# Summary stats for CCS2
## Condition counts
disease_counts <- ccs_validation[,lapply(.SD,sum),.SDcols=names_ccs2]

## AF counts
setDF(ccs_validation)
af_cts <- rep(NA,length(names_ccs2))

for (i in 1:length(names_ccs2)){
  af_cts[i] <- sum(ccs_validation[ccs_validation[,names_ccs2[i]]==1,]$af_5y_sal)
}

## Concordance
# Limit to classes where there are at least 10 AF cases
names_ccs2 <- names_ccs2[af_cts >= 10]
disease_counts <- unlist(disease_counts)[af_cts >= 10]
af_cts <- af_cts[af_cts >= 10]

c <- rep(NA,length(names_ccs2))

for (i in 1:length(names_ccs2)){
  c[i] <- summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=ccs_validation[ccs_validation[,names_ccs2[i]]==1,]))$concordance[1]
}

## Concordance SE (requires bootstrap)
# Bootstrap function
boot <- function(time,status,response,data,runs,size=nrow(data)){
  out <- rep(NA,times=runs)
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE),]
    out[i] <- summary(coxph(Surv(sample[,time],sample[,status]) ~ sample[,response],data=sample))$concordance[1]
  }
  return(out)}
  
# Error robust bootstrap function
boot_err <- function(time,status,response,data,runs,size=nrow(data)){
  out <- rep(NA,times=runs)
  for (i in 1:runs){
    while (is.na(out[i])){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE),]
    mod <- try({coxph(Surv(sample[,time],sample[,status]) ~ sample[,response],data=sample)})
    if (inherits(mod, "try-error")){out[i] <- NA
    } else { out[i] <- summary(mod)$concordance[1] }}
  }
  return(out)
}

# Loop the bootstrap function over each CCS class
se <- rep(NA,length(names_ccs2))
for (i in 1:length(names_ccs2)){
  out <- boot_err(time='af_5y_sal.t',status='af_5y_sal',response='score',
               data=ccs_validation[ccs_validation[,names_ccs2[i]]==1,],runs=200)
  se[i] <- sd(out)
  print(paste0('I just finished subgroup ',i,' out of ',length(names_ccs2),'!'))
}

## Condense output
output <- data.frame(ccs2=names_ccs2,n_disease=unlist(disease_counts),incd_af=af_cts,concordance=c,
                     concordance_se=se)

output$concordance_lb <- output$concordance-1.96*output$concordance_se
output$concordance_ub <- output$concordance+1.96*output$concordance_se

## Save out
write.csv(output,file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/ccs2_output.csv',row.names=F)

###################################################### STEP 5: Plotting
# Load output file with hand-curated disease names (OPTIONAL)
output <- read.csv('/data/arrhythmia/skhurshid/broad_ibm_afrisk/ccs2_output.csv')

# Sort data by concordance (low to high)
output <- output[order(output$concordance),]

# Let's plot the top 50 concordances
png('/data/arrhythmia/skhurshid/broad_ibm_afrisk/concordance_top_ccs2.png',height=700,width=1000)
par(oma=c(1,1,1,1))
par(mar=c(3,30,1,1))
plot(x=output$concordance[64:113],xlim=c(0.700,0.975),y=seq(1,50,1),xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col='#6baed6',cex=1)
axis(1,at=seq(0.700,0.975,0.025),cex=1.4)
axis(2,at=1:50,labels=FALSE,cex=1.4)
mtext('CCS class (level 2)',side=2,line=29.5,cex=1.3)
mtext('Concordance',side=1,line=2.5,cex=1.3)
segments(output$concordance_lb[64:113],1:50,output$concordance_ub[64:113],1:50,col='#6baed6',lwd=2.2)

plot_names <- output$ccs2[64:113]
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i],' (N=',output$n_disease[i],'; AF=',output$incd_af[i],')'),
        side=2,line=1,las=2,at=i)
}

  text(output$concordance_ub[64:113]+0.03,1:50,as.character(round(output$concordance[64:113],3)))
dev.off()

# Let's plot the bottom 50 concordances
png('/data/arrhythmia/skhurshid/broad_ibm_afrisk/concordance_bottom_ccs2.png',height=700,width=1000)
par(oma=c(1,1,1,1))
par(mar=c(3,30,1,1))
plot(x=output$concordance[1:50],xlim=c(0.500,0.800),y=seq(1,50,1),xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col='#e31a1c',cex=1)
axis(1,at=seq(0.500,0.800,0.050),cex=1.4)
axis(2,at=1:50,labels=FALSE,cex=1.4)
mtext('CCS class (level 2)',side=2,line=29.5,cex=1.3)
mtext('Concordance',side=1,line=2.5,cex=1.3)
segments(output$concordance_lb[1:50],1:50,output$concordance_ub[1:50],1:50,col='#e31a1c',lwd=2.2)

plot_names <- output$ccs2[1:50]
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i],' (N=',output$n_disease[i],'; AF=',output$incd_af[i],')'),
        side=2,line=1,las=2,at=i)
}

text(output$concordance_lb[1:50]-0.03,1:50,as.character(round(output$concordance[1:50],3)))
text(0.525,2,as.character(round(output$concordance[2],3)))

dev.off()
############################################################################################################
# CCS LEVEL 3 NOW
############################################################################################################

###################################################### STEP 3: Create prevalent variables based on start_fu
# Gather ccs3 related variables and park in list of names
names_ccs3 <- names(ccs_validation)[c(str_detect(names(ccs_validation),'ccs3') & 
                                        str_detect(names(ccs_validation),'_fin$'))]

str_sub(names_ccs3,-4,-1) <- '_baseline'

# Loop over the list to create baseline variables
## Set DF
setDF(ccs_validation)

## Create list of names for input
input_names <- names_ccs3
str_sub(input_names,-9,-1) <- ''

## Loop over the list (baseline = 1 if CCS condition occurred prior to start_fu)
for (i in 1:length(names_ccs3)){
  ccs_validation[,names_ccs3[i]] <- ifelse(c((ccs_validation[,paste0(input_names[i],'_fin')] == 1) &
                                               (ccs_validation[,paste0(input_names[i],'_fin_d')] <= ccs_validation[,'start_fu'])),1,0)
}

## Set back to DT
setDT(ccs_validation)

###################################################### STEP 4: Concordance analyses
# Summary stats for ccs3
## Condition counts
disease_counts <- ccs_validation[,lapply(.SD,sum),.SDcols=names_ccs3]

## AF counts
setDF(ccs_validation)
af_cts <- rep(NA,length(names_ccs3))

for (i in 1:length(names_ccs3)){
  af_cts[i] <- sum(ccs_validation[ccs_validation[,names_ccs3[i]]==1,]$af_5y_sal)
}

## Concordance
# Limit to classes where there are at least 10 AF cases
names_ccs3 <- names_ccs3[af_cts >= 10]
disease_counts <- unlist(disease_counts)[af_cts >= 10]
af_cts <- af_cts[af_cts >= 10]

c <- rep(NA,length(names_ccs3))

for (i in 1:length(names_ccs3)){
  c[i] <- summary(coxph(Surv(af_5y_sal.t,af_5y_sal) ~ score,data=ccs_validation[ccs_validation[,names_ccs3[i]]==1,]))$concordance[1]
}

## Concordance SE (requires bootstrap)
# Bootstrap function
boot <- function(time,status,response,data,runs,size=nrow(data)){
  out <- rep(NA,times=runs)
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE),]
    out[i] <- summary(coxph(Surv(sample[,time],sample[,status]) ~ sample[,response],data=sample))$concordance[1]
  }
  return(out)}

# Error robust bootstrap function
boot_err <- function(time,status,response,data,runs,size=nrow(data)){
  out <- rep(NA,times=runs)
  for (i in 1:runs){
    while (is.na(out[i])){
      sample <- data[sample(1:nrow(data),size=size,replace=TRUE),]
      mod <- try({coxph(Surv(sample[,time],sample[,status]) ~ sample[,response],data=sample)})
      if (inherits(mod, "try-error")){out[i] <- NA
      } else { out[i] <- summary(mod)$concordance[1] }}
  }
  return(out)
}

# Loop the bootstrap function over each CCS class
se <- rep(NA,length(names_ccs3))
for (i in 1:length(names_ccs3)){
  out <- boot_err(time='af_5y_sal.t',status='af_5y_sal',response='score',
                  data=ccs_validation[ccs_validation[,names_ccs3[i]]==1,],runs=200)
  se[i] <- sd(out)
  print(paste0('I just finished subgroup ',i,' out of ',length(names_ccs3),'!'))
}

## Condense output
output <- data.frame(ccs3=names_ccs3,n_disease=unlist(disease_counts),incd_af=af_cts,concordance=c,
                     concordance_se=se)

output$concordance_lb <- output$concordance-1.96*output$concordance_se
output$concordance_ub <- output$concordance+1.96*output$concordance_se

## Save out
write.csv(output,file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/ccs3_output.csv',row.names=F)

###################################################### STEP 5: Plotting
# Load output file with hand-curated disease names (OPTIONAL)
output <- read.csv('/data/arrhythmia/skhurshid/broad_ibm_afrisk/ccs3_output.csv')

# Sort data by concordance (low to high)
output <- output[order(output$concordance),]

# Let's plot the top 50 concordances
png('/data/arrhythmia/skhurshid/broad_ibm_afrisk/concordance_top_ccs3.png',height=700,width=1000)
par(oma=c(1,1,1,1))
par(mar=c(3,30,1,1))
plot(x=output$concordance[168:217],xlim=c(0.700,0.975),y=seq(1,50,1),xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col='#6baed6',cex=1)
axis(1,at=seq(0.700,0.975,0.025),cex=1.4)
axis(2,at=1:50,labels=FALSE,cex=1.4)
mtext('CCS class (level 3)',side=2,line=29.5,cex=1.3)
mtext('Concordance',side=1,line=2.5,cex=1.3)
segments(output$concordance_lb[168:217],1:50,output$concordance_ub[168:217],1:50,col='#6baed6',lwd=2.2)

plot_names <- output$ccs3[168:217]
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i],' (N=',output$n_disease[i],'; AF=',output$incd_af[i],')'),
        side=2,line=1,las=2,at=i)
}

text(output$concordance_ub[168:217]+0.03,1:50,as.character(round(output$concordance[168:217],3)))
text(0.525,2,as.character(round(output$concordance[2],3)))
dev.off()

# Let's plot the bottom 50 concordances
png('/data/arrhythmia/skhurshid/broad_ibm_afrisk/concordance_bottom_ccs3.png',height=700,width=1000)
par(oma=c(1,1,1,1))
par(mar=c(3,30,1,1))
plot(x=output$concordance[1:50],xlim=c(0.500,0.800),y=seq(1,50,1),xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col='#e31a1c',cex=1)
axis(1,at=seq(0.500,0.800,0.050),cex=1.4)
axis(2,at=1:50,labels=FALSE,cex=1.4)
mtext('CCS class (level 3)',side=2,line=29.5,cex=1.3)
mtext('Concordance',side=1,line=2.5,cex=1.3)
segments(output$concordance_lb[1:50],1:50,output$concordance_ub[1:50],1:50,col='#e31a1c',lwd=2.2)

plot_names <- output$ccs3[1:50]
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i],' (N=',output$n_disease[i],'; AF=',output$incd_af[i],')'),
        side=2,line=1,las=2,at=i)
}

text(output$concordance_lb[c(3:4,7:50)]-0.03,c(3:4,7:50),as.character(round(output$concordance[c(3:4,7:50)],3)))
text(output$concordance_ub[c(1,2,5,6)]+0.03,c(1,2,5,6),as.character(round(output$concordance[c(1,2,5,6)],3)))
text(0.525,18,as.character(round(output$concordance[18],3)))

dev.off()
