# Script to validate Explorys findings in PHS legacy cohort

# Dependencies
library(data.table)
library(stringr)

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

# Load CCS data
load(file='/data/arrhythmia/skhurshid/ccs2/file_2001.RData')

# Load legacy validation
load(file='/data/arrhythmia/skhurshid/ehr_af/vs_021219.RData'); setDT(vs)

# Reduce each only to necessary variables to reduce overhead
ds <- ds[,c('EMPI','start_fu','af_5y_sal','af_5y_sal.t','score','pred5')]
vs <- vs[,c('EMPI','start_fu','af_5y_sal','af_5y_sal.t','score','pred5')]

# Simple bind
all <- rbind(ds,vs)