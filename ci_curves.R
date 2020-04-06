# Script to plot CI curves from Explorys output

# Dependencies
library(data.table)
library(prodlim)
library(survival)
library(Cairo)

# Load data
load("~/Documents/MGH Research/explorys/ci.af_v_pred_risk_EHR_AF_Final.Rdata")

# Plot
CairoPDF(file='~/Documents/MGH Research/explorys/ci_ehr.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(ci.af_v,"cuminc",ylim=c(0,0.15),
     axis2.at=seq(0,0.15,0.05),axis2.las=2,lwd=1.5,background=F,
     atrisk.times=c(0,1,2,3,4,5),col=c("#800026","#fd8d3c",'#fed976'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.15,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("5-yr risk \u2265 5%","5-yr risk 2.5-5%","5-yr risk <2.5%"),
     atrisk.title=("                     "),atrisk.pos=-0.35,atrisk.line=c(1.2,2.5,3.8),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative incidence of AF (%)",side=2,line=-2.5,at=0.075,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('AF risk level',side=1, line=-0.5,cex=1.8,at=-0.35)
dev.off()