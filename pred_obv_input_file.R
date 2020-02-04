# Script to make pred/obv plot with Explorys input data
# Dependencies
library(data.table)

# Load input data
output <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/pred_obv_input.csv')

# Plotting
# Open file connection
png("/data/arrhythmia/skhurshid/broad_ibm_afrisk/pred_obv.png",height=960,width=1800,res=100)

## Params
par(oma=c(1,1,0,0))

## x-scales
x1 <- output[name=='CHADS-VASc']$x1
x2 <- output[name=='CHARGE-AF']$x1
x3 <- output[name=='EHR-AF']$x1
x4 <- output[name=='C2HEST']$x1

## y-values
chads <- output[name=='CHADS-VASc']$x2
charge <- output[name=='CHARGE-AF']$x2
ehraf <- output[name=='EHR-AF']$x2
chest <- output[name=='C2HEST']$x2

## colors
col1 <- "#1f78b4"; col2 <- "#33a02c"; col3 <- "#e31a1c"; col4 <- "#ff7f00"

# Plot
## Point estimates
plot(x=c(x1,x2,x3,x4),y=c(chads,charge,ehraf,chest),
     pch=19,col=c(rep(col1,length(x1)),rep(col2,length(x2)),rep(col3,length(x3)),rep(col4,length(x4))),
     xlim=c(-0.5,26),ylim=c(0,30),
     xaxt='n',yaxt='n',xlab='',ylab='',cex=1.6,frame.plot = F)

## CIs
#segments(x3,ehraf_obv[,2],x3,ehraf_obv[,3],col=col3)
#segments(x2,charge_obv[,2],x2,charge_obv[,3],col=col2)
#segments(x1,chads_obv[,2],x1,chads_obv[,3],col=col1)
#segments(x4,chest_obv[,2],x4,chest_obv[,3],col=col4)

## Axes
axis(side=1,cex.axis=1.5,at=seq(0,25,1),labels=seq(0,25,1),pos=-0.5)
axis(side=2,cex.axis=1.5,at=seq(0,30,5),labels=seq(0,30,5),las=1,pos=-0.5)

## Labels
mtext(side=2,"Observed 5-Year AF (%)",line=1.8,cex=1.6)
mtext(side=1,"Predicted 5-Year AF Risk (%)",line=3,cex=1.6)

## Segments
segments(0,-0.5,-0.5,-0.5)
segments(-0.5,0,-0.5,-0.5)
segments(25,-0.5,25.5,-0.5)

# 45 Line
black_trans <- adjustcolor('black', alpha.f = 0.6) 
segments(-0.5,-0.5,25,25,col=black_trans)

## Legend
legend(-0.2,29.5,c('CHADS-VASc','CHARGE-AF','EHR-AF','C2HEST'),
       pch=19,col=c(col1,col2,col3,col4),cex=1.5,bty='n')

dev.off()
