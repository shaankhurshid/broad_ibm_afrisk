# Script to make pred/obv plot with Explorys input data
# Dependencies
library(data.table)

# Load input data
coord <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/pred_obv_input2.csv')
input <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/pred_obv_input.csv')

# Plotting
# Open file connection
pdf("/data/arrhythmia/skhurshid/broad_ibm_afrisk/pred_obv_lines.pdf",height=5,width=5,
    pointsize=1)

## Params
par(oma=c(1,1,0,0))

## x-scales
x1 <- coord[V3=='CHADS-VASc']$x1
x2 <- coord[V3=='CHARGE-AF']$x1
x3 <- coord[V3=='EHR-AF']$x1
x4 <- coord[V3=='C2HEST']$x1

## y-values
chads <- input[V4=='CHADS-VASc']$V1
charge <- input[V4=='CHARGE-AF']$V1
ehraf <- input[V4=='EHR-AF']$V1
chest <- input[V4=='C2HEST']$V1

## colors
col1 <- "#1f78b4"; col2 <- "#33a02c"; col3 <- "#e31a1c"; col4 <- "#ff7f00"

# Plot
## Point estimates
plot(x=c(x1,x2,x3,x4),y=c(chads,charge,ehraf,chest),
     pch=19,
     col=c(rep(col1,length(x1)),rep(col2,length(x2)),rep(col3,length(x3)),rep(col4,length(x4))),
     xlim=c(-0.5,26),ylim=c(0,25),
     xaxt='n',yaxt='n',xlab='',ylab='',cex=1.4,frame.plot = F)

## CIs
segments(x1,input[V4=='CHADS-VASc']$V2,x1,input[V4=='CHADS-VASc']$V3,col=col1)
segments(x2,input[V4=='CHARGE-AF']$V2,x2,input[V4=='CHARGE-AF']$V3,col=col2)
segments(x3,input[V4=='EHR-AF']$V2,x3,input[V4=='EHR-AF']$V3,col=col3)
segments(x4,input[V4=='C2HEST']$V2,x4,input[V4=='C2HEST']$V3,col=col4)

## Axes
axis(side=1,cex.axis=1.5,at=seq(0,25,5),labels=seq(0,25,5),pos=-0.5)
axis(side=2,cex.axis=1.5,at=seq(0,25,5),labels=seq(0,25,5),las=1,pos=-0.5)

## Labels
mtext(side=2,"Observed 5-Year AF (%)",line=1.8,cex=1.6)
mtext(side=1,"Predicted 5-Year AF Risk (%)",line=3,cex=1.6)

## Segments
segments(0,-0.5,-0.5,-0.5)
segments(-0.5,0,-0.5,-0.5)
segments(25,-0.5,25.5,-0.5)

## Connecting lines
for (i in 1:(length(x1)-1)){
  segments(x1[i],chads[i],x1[i+1],chads[i+1],col=col1)
}

for (i in 1:(length(x2)-1)){
segments(x2[i],charge[i],x2[i+1],charge[i+1],col=col2)
}

for (i in 1:(length(x3)-1)){
segments(x3[i],ehraf[i],x3[i+1],ehraf[i+1],col=col3)
}

for (i in 1:(length(x4)-1)){
segments(x4[i],chest[i],x4[i+1],chest[i+1],col=col4)
}

# 45 Line
black_trans <- adjustcolor('black', alpha.f = 0.6) 
segments(-0.5,-0.5,25,25,col=black_trans,lwd=2,lty=2)

## Legend
legend(-0.2,24.5,c('CHADS-VASc','CHARGE-AF','EHR-AF','C2HEST'),
       pch=19,col=c(col1,col2,col3,col4),cex=1.5,bty='n')

dev.off()
