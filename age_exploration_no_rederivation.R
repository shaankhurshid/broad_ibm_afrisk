# Dependencies
library(data.table); library(survival)

# Load data
load(file='/data/arrhythmia/skhurshid/ehr_af/vs_032120.RData')

# Script to explore coefficients in subgroups

# Define explore function
explore_age <- function(time,status,age_variable,min_age,max_age,
                        age_step,data,risk_score,risk,threshold){
i <- 1
out <- list()
for (age in seq(min_age,max_age,age_step)){
  if(age == max_age-age_step){
  subset <- data[c(data[,age_variable] >= age & data[,age_variable] <= age+age_step),]
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
  out[[i]] <- as.numeric(c(paste0(age),n_af,n_total,ci,
                           as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                           as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3]),
                           ppv,npv,sens,spec))
  names(out[[i]]) <- c('age','n_af','n_total','ci','ci_lower','ci_upper',
                       'c_stat','c_stat_lb','c_stat_ub','cal','cal_lb','cal_ub',
                       'ppv','ppv_lower','ppv_upper',
                       'npv','npv_lower','npv_upper',
                       'sens','sens_lower','sens_upper',
                       'spec','spec_lower','spec_upper')
  print(paste0('Just finished model ',i,' out of ',((max_age-min_age)/age_step),'!'))
  break
  } else {
    subset <- data[c(data[,age_variable] >= age & data[,age_variable] < age+age_step),]
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
    out[[i]] <- as.numeric(c(paste0(age),n_af,n_total,ci,
                             as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                             as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3]),
                             ppv,npv,sens,spec))
    names(out[[i]]) <- c('age','n_af','n_total','ci','ci_lower','ci_upper',
                         'c_stat','c_stat_lb','c_stat_ub','cal','cal_lb','cal_ub',
                         'ppv','ppv_lower','ppv_upper',
                         'npv','npv_lower','npv_upper',
                         'sens','sens_lower','sens_upper',
                         'spec','spec_lower','spec_upper')
    print(paste0('Just finished model ',i,' out of ',((max_age-min_age)/age_step),'!'))
    i <- i+1
  }
}
return(data.frame(do.call(rbind,out)))
}

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

# Run explore functions
output_age <- explore_age(time='af_5y_sal.t',status='af_5y_sal',min_age=45,max_age=95,age_step=1,
                      age_variable='start_fu_age',data=vs,risk_score='score',risk='pred5',threshold=0.05)

output_sex <- explore_categorical(time='af_5y_sal.t',status='af_5y_sal',variable='Gender',
                              data=vs,risk_score='score',risk='pred5',threshold=0.05)

# Discrimination
pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/cstat_age.pdf',height=3,width=9,
    pointsize=5)
par(mar=c(3,2,1,1),oma=c(3,2,1,1))
plot(x=output$age,y=output$c_stat,col='#f46d43',xaxt='n',yaxt='n',xlab='',ylab='',
     pch=19,cex=1,ylim=c(0.3,1),xlim=c(45,95),bty='n')

axis(1,at=seq(45,95,5),cex.axis=1.2)
mtext("Age",1,line=3,cex=1.5)

axis(2,at=seq(0.3,1,0.1),las=2,cex.axis=1.2,pos=44.5)
mtext("C-Statistic",2,line=2,cex=1.5)

text(x=output$age,y=output$c_stat_ub+0.02,paste0(output$n_af))

segments(x0=output$age,y0=output$c_stat_lb,x1=output$age,y1=output$c_stat_ub,col='#f46d43')
dev.off()

# Calibration
pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/cal_age.pdf',height=3,width=9,
    pointsize=5)
par(mar=c(3,2,1,1),oma=c(3,2,1,1))
plot(x=output$age,y=output$cal,col='#1a9850',xaxt='n',yaxt='n',xlab='',ylab='',
     pch=19,cex=1,ylim=c(-2,2.5),xlim=c(45,95),bty='n')

axis(1,at=seq(45,95,5),cex.axis=1.2)
mtext("Age",1,line=3,cex=1.5)

axis(2,at=seq(-2,2.5,0.5),las=2,pos=44.5,cex.axis=1.2)
mtext("Calibration Slope",2,line=2,cex=1.5)

text(x=output$age,y=output$cal_ub+0.2,paste0(output$n_af))

segments(x0=output$age,y0=output$cal_lb,x1=output$age,y1=output$cal_ub,col='#1a9850')
segments(45,1,95,1,lty=5,col='black')
dev.off()

# CI
pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/ci_age.pdf',height=3,width=9,
    pointsize=5)
par(mar=c(3,2,1,1),oma=c(3,2,1,1))
plot(x=output$age,y=output$ci,col='#2b8cbe',xaxt='n',yaxt='n',xlab='',ylab='',
     pch=19,cex=1,ylim=c(0,40),xlim=c(45,95),bty='n')

axis(1,at=seq(45,95,5),cex.axis=1.2)
mtext("Age",1,line=3,cex=1.5)

axis(2,at=seq(0,40,5),las=2,pos=44.5,cex.axis=1.2)
mtext("Cumulative Incidence of AF (%)",2,line=2,cex=1.5)

text(x=output$age,y=output$ci_upper+1.5,paste0(output$n_af))

segments(x0=output$age,y0=output$ci_lower,x1=output$age,y1=output$ci_upper,col='#2b8cbe')
dev.off()

# Sensitivity
pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/sens_age.pdf',height=3,width=9,
    pointsize=5)
par(mar=c(3,2,1,1),oma=c(3,2,1,1))
plot(x=output$age,y=output$sens,col='#f46d43',xaxt='n',yaxt='n',xlab='',ylab='',
     pch=19,cex=1,ylim=c(0,1),xlim=c(45,95),bty='n')

axis(1,at=seq(45,95,5),cex.axis=1.2)
mtext("Age",1,line=3,cex=1.5)

axis(2,at=seq(0,1,0.1),las=2,cex.axis=1.2,pos=44.5)
mtext("Sensitivity for AF at 5% risk threshold (%)",2,line=2,cex=1.5)

par(xpd=TRUE)
for (i in 1:length(output$age)){
  if (is.na(output$sens_upper[i])){
    text(x=output$age[i],y=output$sens[i]+0.04,paste0(output$n_af[i]))
  } else
    text(x=output$age[i],y=output$sens_upper[i]+0.06,paste0(output$n_af[i]))
}

segments(x0=output$age,y0=output$sens_lower,x1=output$age,y1=output$sens_upper,col='#f46d43')
dev.off()

# Specificity
pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/spec_age.pdf',height=3,width=9,
    pointsize=5)
par(mar=c(3,2,1,1),oma=c(3,2,1,1))
plot(x=output$age,y=output$spec,col='#1a9850',xaxt='n',yaxt='n',xlab='',ylab='',
     pch=19,cex=1,ylim=c(0,1),xlim=c(45,95),bty='n')

axis(1,at=seq(45,95,5),cex.axis=1.2)
mtext("Age",1,line=3,cex=1.5)

axis(2,at=seq(0,1,0.1),las=2,cex.axis=1.2,pos=44.5)
mtext("Specificity for AF at 5% risk threshold (%)",2,line=2,cex=1.5)

par(xpd=TRUE)
for (i in 1:length(output$age)){
  if (is.na(output$spec_upper[i])){
    text(x=output$age[i],y=output$spec[i]+0.04,paste0(output$n_af[i]))
  } else
    text(x=output$age[i],y=output$spec_upper[i]+0.06,paste0(output$n_af[i]))
}

segments(x0=output$age,y0=output$spec_lower,x1=output$age,y1=output$spec_upper,col='#1a9850')
dev.off()