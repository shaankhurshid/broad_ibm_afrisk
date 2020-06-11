# Script to produce NRI

# Dependencies
library(data.table)

# Define binary cutoff
mydata[,':='(score_pred5_above5 = ifelse(score_pred5 >= 5,1,0))] # Repeat for each score using predicted risk variable


# Replace time with time variable, event with AF variable
# t0 should be 5 years (change to days if that is how you are coding it)

# NRIs
ehr_v_charge <- nricens(p.std=mydata$charge_pred5_above5, p.new=mydata$predict_af_pred5_above5,
                              time=mydata$incd_af_5y.t, event=mydata$incd_af_5y, cut=0.5,
                              niter = 10, t0=5,updown='category')

ehr_v_chest <- nricens(p.std=mydata$chest_pred5_above5, p.new=mydata$predict_af_pred5_above5,
                        time=mydata$incd_af_5y.t, event=mydata$incd_af_5y, cut=0.5,
                        niter = 10, t0=5,updown='category')

ehr_v_chads <- nricens(p.std=mydata$chads_pred5_above5, p.new=mydata$predict_af_pred5_above5,
                        time=mydata$incd_af_5y.t, event=mydata$incd_af_5y, cut=0.5,
                        niter = 10, t0=5,updown='category')