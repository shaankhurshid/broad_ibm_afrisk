# Depends
library(data.table)
library(nricens)

# Load original and recalibrated datasets
charge_age_discrim <- fread(file='/data/arrhythmia/skhurshid/heterogeneity/charge_output_age.csv')
charge_age_discrim_recal <- fread(file='/data/arrhythmia/skhurshid/heterogeneity/charge_output_age_recal.csv')

# Function to generate predicted risks based on a subset
calibrated_risk <- function(subset,variables_list,time,status,output_name='young_pred5'){
  new_formula <- formula(paste0('Surv(subset[,get(time)],subset[,get(status)]) ~ ',paste0(variables_list,collapse=' + ')))
  mod <- coxph(new_formula,data=subset)
  subset$lp <- predict(mod,type='lp')
  avg_beta <- mean(subset$lp)
  res <- coxph(Surv(subset[,get(time)],subset[,get(status)]) ~ lp, data=subset)
  km <- survfit(res, data=data.frame(x1=mean(lp)),type="kaplan-meier")
  s0 <- summary(km, times=c(5))$surv
  subset[,paste0(output_name) := (1-(s0)^exp(lp - (avg_beta)))]
}

charge_vars <- c('start_fu_age_5','race_binary','ht_cm_atStartFu_10','wt_kg_atStartFu_15','sbp_atStartFu_20',
                 'dbp_atStartFu_10','tobacco_fin_prevAtstartFu','htn_fin_prevAtstartFu','dm_fin_prevAtstartFu',
                 'heartFailure_fin_prevAtstartFu','mi_fin_prevAtstartFu')

### CASE 1: Young people
# Define subset
young <- vs[c(start_fu_age >= 45 & start_fu_age < 50)]
young <- calibrated_risk(subset=young,variables_list=charge_vars,time='af_5y_sal.t',status='af_5y_sal',output_name='young_pred5')

# NRI
young_nri <- nricens(p.std=young$charge.pred5, p.new=young$young_pred5,
                        time=young$af_5y_sal.t, event=young$af_5y_sal, cut=0.025,
                        niter = 10, t0=5,updown='category')

### CASE 2: Old people
# Define subset
old <- vs[c(start_fu_age >= 80 & start_fu_age < 85)]
old <- calibrated_risk(subset=old,variables_list=charge_vars,time='af_5y_sal.t',status='af_5y_sal',output_name='old_pred5')

# NRI
old_nri <- nricens(p.std=old$charge.pred5, p.new=old$old_pred5,
                 time=old$af_5y_sal.t, event=old$af_5y_sal, cut=0.25,
                 niter = 10, t0=5,updown='category')


### CASE 3: HF
# Define subset
hf <- vs[heartFailure_fin_prevAtstartFu==1]
hf <- calibrated_risk(subset=hf,variables_list=charge_vars,time='af_5y_sal.t',status='af_5y_sal',output_name='hf_pred5')

# NRI
hf_nri <- nricens(p.std=hf$charge.pred5, p.new=hf$hf_pred5,
                   time=hf$af_5y_sal.t, event=hf$af_5y_sal, cut=0.15,
                   niter = 10, t0=5,updown='category')