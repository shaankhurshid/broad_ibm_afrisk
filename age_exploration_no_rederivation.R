# Dependencies
library(data.table); library(survival)

# Script to explore coefficients in subgroups

# Define explore function
explore_age <- function(time,status,age_variable,min_age,max_age,age_step,data,risk_score){
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
  out[[i]] <- data.frame(matrix(ncol=12,nrow=0))
  out[[i]] <- as.numeric(c(paste0(age),n_af,n_total,ci,
                           as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                           as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3])))
  names(out[[i]]) <- c('age','n_af','n_total','ci','ci_lower','ci_upper',
                       'c_stat','c_stat_lb','c_stat_ub','cal','cal_lb','cal_ub')
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
    out[[i]] <- data.frame(matrix(ncol=12,nrow=0))
    out[[i]] <- as.numeric(c(paste0(age),n_af,n_total,ci,
                  as.numeric(summary(mod)$concordance[1]),as.numeric(summary(mod)$concordance[1]-1.96*summary(mod)$concordance[2]),as.numeric(summary(mod)$concordance[1]+1.96*summary(mod)$concordance[2]),
                  as.numeric(mod$coefficients[1]),as.numeric(mod$coefficients[1]-1.96*summary(mod)$coefficients[3]),as.numeric(mod$coefficients[1]+1.96*summary(mod)$coefficients[3])))
    names(out[[i]]) <- c('age','n_af','n_total','ci','ci_lower','ci_upper',
                         'c_stat','c_stat_lb','c_stat_ub','cal','cal_lb','cal_ub')
    print(paste0('Just finished model ',i,' out of ',((max_age-min_age)/age_step),'!'))
    i <- i+1
  }
}
return(do.call(rbind,out))
}

# Run explore function
data <- vs # put your time variable in time, AF variable in status, baseline age variable in age_variable

output <- data.frame(explore_age(time='af_5y_sal.t',status='af_5y_sal',min_age=45,max_age=95,age_step=1,
                      age_variable='start_fu_age',data=data,risk_score='score'))

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