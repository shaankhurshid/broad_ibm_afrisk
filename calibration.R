# Script for generating calibration plots with predicted risk distribution on top

# Dependencies
library(rms)

## CHARGE
survmod<-with(data,Surv(af_status,af_time))
fit<-cph(survmod~charge,data=data,surv=TRUE,time.inc=5*365.25,u=5*365.25,x=T,y=T)
fit2 <- coxph(survmod~charge,data=data)
cal_charge <- c(fit2$coefficients[1],confint(fit2,"charge")[1],confint(fit2,"charge")[2]) # beta and 95%CI corresponds to calibration slope

# plot calibrations
cal.charge<-calibrate(fit,u=5*365.25,cmethod='hare',B=10) # ultimately want B=200 but takes a while, B=10 good enough to see how it looks

pdf('charge.pdf',height=3,width=3,pointsize=3)
par(oma=c(1,1,1,1))
par(oma=c(1,1,1,1))

col=paste(rgb(252,146,114,maxColorValue=255),sep="")

plot(cal.charge,scat1d.opts=list(frac=0.1,side=1),xlim=c(1,0.2),ylim=c(1,0.2),xaxt="n",yaxt="n", xlab="Predicted 5-year risk of AF", ylab="Proportion with AF at 5 years", subtitles=F, par.corrected=list(col=col2, lty=1, lwd=2), bty='n', cex.lab=1.25)
axis(1,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))
axis(2,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0), las=1)
mtext(text="CHARGE-AF", side=3, line=-0.5, adj=0.5, cex=1.2)
mtext(text="Calibration slope", side=3, line=-2.5, at=0.85, cex=0.9)
mtext(text="0.99 (0.96-1.01)", side=3, line=-3.5, at=0.85, cex=0.9)

legend(0.5,0.8,c('optimal','observed','optimism-corrected'),col=c('darkgray','black',col),lty=1,lwd=1.5,bty='n',cex=0.75)
segments(0,0,0.9,0.9,col='darkgray',lty=1)

dev.off()

## ehraf
survmod<-with(data,Surv(af_status,af_time))
fit<-cph(survmod~ehraf,data=data,surv=TRUE,time.inc=5*365.25,u=5*365.25,x=T,y=T)
fit2 <- coxph(survmod~ehraf,data=data)
cal_ehraf <- c(fit2$coefficients[1],confint(fit2,"ehraf")[1],confint(fit2,"ehraf")[2]) # beta and 95%CI corresponds to calibration slope

# plot calibrations
cal.ehraf<-calibrate(fit,u=5*365.25,cmethod='hare',B=10) # ultimately want B=200 but takes a while, B=10 good enough to see how it looks

pdf('ehraf.pdf',height=3,width=3,pointsize=3)
par(oma=c(1,1,1,1))
par(oma=c(1,1,1,1))

plot(cal.ehraf,scat1d.opts=list(frac=0.1,side=1),xlim=c(1,0.2),ylim=c(1,0.2),xaxt="n",yaxt="n", xlab="Predicted 5-year risk of AF", ylab="Proportion with AF at 5 years", subtitles=F, par.corrected=list(col=col2, lty=1, lwd=2), bty='n', cex.lab=1.25)
axis(1,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))
axis(2,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0), las=1)
mtext(text="EHR-AF", side=3, line=-0.5, adj=0.5, cex=1.2)
mtext(text="Calibration slope", side=3, line=-2.5, at=0.85, cex=0.9)
mtext(text="0.99 (0.96-1.01)", side=3, line=-3.5, at=0.85, cex=0.9)

legend(0.5,0.8,c('optimal','observed','optimism-corrected'),col=c('darkgray','black',col),lty=1,lwd=1.5,bty='n',cex=0.75)
segments(0,0,0.9,0.9,col='darkgray',lty=1)

dev.off()

## CHADS
survmod<-with(data,Surv(af_status,af_time))
fit<-cph(survmod~chads,data=data,surv=TRUE,time.inc=5*365.25,u=5*365.25,x=T,y=T)
fit2 <- coxph(survmod~chads,data=data)
cal_chads <- c(fit2$coefficients[1],confint(fit2,"chads")[1],confint(fit2,"chads")[2]) # beta and 95%CI corresponds to calibration slope

# plot calibrations
cal.chads<-calibrate(fit,u=5*365.25,cmethod='hare',B=10) # ultimately want B=200 but takes a while, B=10 good enough to see how it looks

pdf('chads.pdf',height=3,width=3,pointsize=3)
par(oma=c(1,1,1,1))
par(oma=c(1,1,1,1))

plot(cal.chads,scat1d.opts=list(frac=0.1,side=1),xlim=c(1,0.2),ylim=c(1,0.2),xaxt="n",yaxt="n", xlab="Predicted 5-year risk of AF", ylab="Proportion with AF at 5 years", subtitles=F, par.corrected=list(col=col2, lty=1, lwd=2), bty='n', cex.lab=1.25)
axis(1,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))
axis(2,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0), las=1)
mtext(text=expression(paste(plain("CHA") [plain("2")], plain("DS") [plain("2")], plain("-VASc"))), side=3, line=-0.5, adj=0.5, cex=1.2)
mtext(text="Calibration slope", side=3, line=-2.5, at=0.85, cex=0.9)
mtext(text="0.99 (0.96-1.01)", side=3, line=-3.5, at=0.85, cex=0.9)

legend(0.5,0.8,c('optimal','observed','optimism-corrected'),col=c('darkgray','black',col),lty=1,lwd=1.5,bty='n',cex=0.75)
segments(0,0,0.9,0.9,col='darkgray',lty=1)

dev.off()

## CHARGE
survmod<-with(data,Surv(af_status,af_time))
fit<-cph(survmod~chest,data=data,surv=TRUE,time.inc=5*365.25,u=5*365.25,x=T,y=T)
fit2 <- coxph(survmod~chest,data=data)
cal_chest <- c(fit2$coefficients[1],confint(fit2,"chest")[1],confint(fit2,"chest")[2]) # beta and 95%CI corresponds to calibration slope

# plot calibrations
cal.chest<-calibrate(fit,u=5*365.25,cmethod='hare',B=10) # ultimately want B=200 but takes a while, B=10 good enough to see how it looks

pdf('chest.pdf',height=3,width=3,pointsize=3)
par(oma=c(1,1,1,1))
par(oma=c(1,1,1,1))

plot(cal.chest,scat1d.opts=list(frac=0.1,side=1),xlim=c(1,0.2),ylim=c(1,0.2),xaxt="n",yaxt="n", xlab="Predicted 5-year risk of AF", ylab="Proportion with AF at 5 years", subtitles=F, par.corrected=list(col=col2, lty=1, lwd=2), bty='n', cex.lab=1.25)
axis(1,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))
axis(2,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0), las=1)
mtext(text=expression(paste(plain("C") [plain("2")], plain("HEST"))), side=3, line=-0.5, adj=0.5, cex=1.2)
mtext(text="Calibration slope", side=3, line=-2.5, at=0.85, cex=0.9)
mtext(text="0.99 (0.96-1.01)", side=3, line=-3.5, at=0.85, cex=0.9)

legend(0.5,0.8,c('optimal','observed','optimism-corrected'),col=c('darkgray','black',col),lty=1,lwd=1.5,bty='n',cex=0.75)
segments(0,0,0.9,0.9,col='darkgray',lty=1)

dev.off()