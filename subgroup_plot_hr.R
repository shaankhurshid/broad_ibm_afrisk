# Dependenices
library(data.table)

# Manual input
f <- c(3.175,3.15,3.201)
m <- c(2.605,2.586,2.624)
hf <- c(1.572, 1.548, 1.596)
stroke <- c(2.142,2.104,2.181)
all <- c(2.87, 2.86, 2.89)

# Bind results
ehr <- data.frame(rbind(all,hf,stroke,m,f))
ehr$subgroup <- c('Overall','Heart failure','Stroke','Males','Females')

# Plot
pdf('/data/arrhythmia/skhurshid/broad_ibm_afrisk/subgroup_hr.pdf',height=1.5,width=4.4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,7.5,1,1))
col <- c('black','#1f78b4','#33a02c','#ff7f00','#e31a1c')
plot(x=ehr[,1],y=seq(1,5,1),xlim=c(1,3.50),ylim=c(0.75,5),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=c(18,rep(19,4)),col=col,cex=c(3,rep(1.4,4)),bty='n')
axis(1,at=seq(1,3.50,0.25),cex=1.4,pos=0)
axis(2,at=1:5,labels=FALSE,cex=1.4)
mtext('Subgroup',side=2,line=7,cex=1.3)
mtext('Standardized Hazard Ratio',side=1,line=3.5,cex=1.3)
segments(ehr$X2,1:5,ehr$X3,1:5,col=col,lwd=2.2)

plot_names <- c('Overall','Heart failure','Stroke','Male','Female')
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i]),
        side=2,line=1,las=2,at=i)
}

dev.off()
