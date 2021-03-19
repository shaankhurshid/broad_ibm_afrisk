# Plotting discrimination versus age for CHARGE-AF

# Depends
library(data.table)

# Load data
charge_age_discrim <- fread(file='/data/arrhythmia/skhurshid/heterogeneity/explorys_charge_output_age.csv')

# Plot discrimination versus age
pdf(file='/data/arrhythmia/skhurshid/heterogeneity/explorys_charge_age_discrim.pdf',height=4,width=8,pointsize=5)
par(mar=c(6,1,1,1),oma=c(1,1,1,1))

y_variable <- charge_age_discrim$c_stat
col <- c(rep('#e31a1c',length(y_variable)-1),'darkgray')
        
plot(x=1:length(y_variable),y=y_variable,col=col,
     bty='n',xaxt='n',yaxt='n',pch=19,xlim=c(0,length(y_variable)),
     ylim=c(0.5,0.8),xlab='',ylab='',cex=2)

axis(1,at=1:length(y_variable),cex.axis=2,
     labels=c(paste0(as.character(as.numeric(charge_age_discrim$age[1:(length(y_variable)-1)])),'-',
                     as.character(as.numeric(charge_age_discrim$age[1:(length(y_variable)-1)])+4)),'All'))
axis(2,cex.axis=2,at=seq(0.5,0.8,0.05),las=2,pos=0.6)

mtext("Concordance index",2,line=-1,cex=2.5)
mtext("Age",1,line=4,cex=2.5)

text(1:length(y_variable),charge_age_discrim$c_stat_ub+0.015,as.character(charge_age_discrim$n_event),cex=1.2)

segments(1:length(y_variable),charge_age_discrim$c_stat_lb,
         1:length(y_variable),charge_age_discrim$c_stat_ub,
         col=col,lwd=1.5)

dev.off()

# Plot calibration slope versus age
pdf(file='/data/arrhythmia/skhurshid/heterogeneity/charge_age_cal_slope.pdf',height=4,width=8,pointsize=5)
par(mar=c(6,1,1,1),oma=c(1,1,1,1))

y_variable <- charge_age_discrim$cal
col <- c(rep('#e34a33',length(y_variable)-1),'darkgray')

plot(x=1:length(y_variable),y=y_variable,col=col,
     bty='n',xaxt='n',yaxt='n',pch=19,xlim=c(0,length(y_variable)),
     ylim=c(0,2),xlab='',ylab='',cex=2)

axis(1,at=1:length(y_variable),cex.axis=2,
     labels=c(paste0(as.character(as.numeric(charge_age_discrim$age[1:(length(y_variable)-1)])),'-',
                     as.character(as.numeric(charge_age_discrim$age[1:(length(y_variable)-1)])+4)),'All'))
axis(2,cex.axis=2,at=seq(0,2,0.2),las=2,pos=0.6)

mtext("Calibration slope",2,line=-1,cex=2.5)
mtext("Age",1,line=4,cex=2.5)

segments(0.6,1,length(y_variable),1,lty=5,col='black')

text(1:length(y_variable),charge_age_discrim$cal_ub+0.1,as.character(charge_age_discrim$n_event),cex=1.2)

segments(1:length(y_variable),charge_age_discrim$cal_lb,
         1:length(y_variable),charge_age_discrim$cal_ub,
         col=col,lwd=1.5)

dev.off()


# Plot relative prediction error versus age
pdf(file='/data/arrhythmia/skhurshid/heterogeneity/charge_age_rel_error.pdf',height=4,width=8,pointsize=5)
par(mar=c(6,1,1,1),oma=c(1,1,1,1))

y_variable <- charge_age_discrim$rel_error_mean*100
col <- c(rep('#31a354',length(y_variable)-1),'darkgray')

plot(x=1:length(y_variable),y=y_variable,col=col,
     bty='n',xaxt='n',yaxt='n',pch=19,xlim=c(0,length(y_variable)),
     ylim=c(0,100),xlab='',ylab='',cex=2)

axis(1,at=1:length(y_variable),cex.axis=2,
     labels=c(paste0(as.character(as.numeric(charge_age_discrim$age[1:(length(y_variable)-1)])),'-',
                     as.character(as.numeric(charge_age_discrim$age[1:(length(y_variable)-1)])+4)),'All'))
axis(2,cex.axis=2,at=seq(0,100,10),las=2,pos=0.6)

mtext("Relative prediction error (%)",2,line=-1,cex=2.5)
mtext("Age",1,line=4,cex=2.5)

text(1:length(y_variable),y_variable+10,as.character(charge_age_discrim$n_event),cex=1.2)

dev.off()