# Dependencies
library(data.table)

# Load decile points from Explorys
ehraf <- fread(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/cal_qual_ehr_af_raw_values.csv')
charge <- fread(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/cal_qual_charge_raw_values.csv')

# Plot 1: CHARGE-AF
pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/charge_uncalibrated.pdf',
    height=4,width=4,pointsize=5)
par(mar=c(3,3,2,1),oma=c(3,3,1,1))
plot(x=charge$x,y=charge$y,xlim=c(0,35),ylim=c(0,35),xaxt='n',yaxt='n',
     xlab='',ylab='',bty='n',col='#e31a1c',pch=19)

axis(1,at=seq(0,35,5),cex.axis=1.4)
axis(2,at=seq(0,35,5),cex.axis=1.4,las=2)

mtext("Predicted 5-year risk of AF",1,cex=1.8,line=3)
mtext("Proportion with AF at 5 years",2,cex=1.8,line=3)
mtext("CHARGE-AF",3,cex=2,line=0)

segments(0,0,36,36,col='darkgray',lty=5)

## Connecting lines
for (i in 1:(length(charge$x)-1)){
  segments(charge$x[i],charge$y[i],charge$x[i+1],charge$y[i+1],col='#e31a1c')
}

dev.off()

# Plot 1: EHR-AF
pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/ehraf_uncalibrated.pdf',
    height=4,width=4,pointsize=5)
par(mar=c(3,3,2,1),oma=c(3,3,1,1))
plot(x=ehraf$x,y=ehraf$y,xlim=c(0,35),ylim=c(0,35),xaxt='n',yaxt='n',
     xlab='',ylab='',bty='n',col='#4575b4',pch=19)

axis(1,at=seq(0,35,5),cex.axis=1.4)
axis(2,at=seq(0,35,5),cex.axis=1.4,las=2)

mtext("Predicted 5-year risk of AF",1,cex=1.8,line=3)
mtext("Proportion with AF at 5 years",2,cex=1.8,line=3)
mtext("EHR-AF",3,cex=2,line=0)

segments(0,0,36,36,col='darkgray',lty=5)

## Connecting lines
for (i in 1:(length(ehraf$x)-1)){
  segments(ehraf$x[i],ehraf$y[i],ehraf$x[i+1],ehraf$y[i+1],col='#4575b4')
}

dev.off()

