# This is a script to plot concordance waterfalls

###################################################### Plotting
# Plot 1: Partners validation set
partners <- read.csv('ccs1_partners.csv') # Input of plot A

# Sort data by concordance (low to high)
partners <- partners[order(partners$concordance),]

png('concordance_ccs1_phs.png',height=500,width=800)
par(oma=c(1,1,1,1))
par(mar=c(3,25,1,1))
col <- ifelse(c((round(mean(partners$concordance),3) >= partners$concordance_lb) & (round(mean(partners$concordance),3) <= partners$concordance_ub)),'black',
           ifelse(partners$concordance_lb > mean(partners$concordance),'#6baed6','#e34a33'))

plot(x=partners$concordance,y=seq(1,18,1),xlim=c(0.700,1.00),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col=col,cex=1.5)

axis(1,at=seq(0.700,1.00,0.050),cex.axis=1.6)
axis(2,at=1:18,labels=FALSE,cex.axis=1.6)
mtext('CCS class (level 1)',side=2,line=24.5,cex=1.5)
mtext('Concordance',side=1,line=2.5,cex=1.5)

segments(partners$concordance_lb,1:18,partners$concordance_ub,1:18,col=col,lwd=2.2)

legend(0.885,4.2,legend=c('Superior','Expected','Inferior'),lty=1,lwd=3,
       col=c('#6baed6','black','#e34a33'),bty='n',cex=1.5)

plot_names <- partners$ccs1
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i],' (N=',partners$n_disease[i],'; AF=',partners$incd_af[i],')'),
        side=2,line=1,las=2,at=i,cex=1.2)
}

text(partners$concordance_ub+0.03,1:17,as.character(round(partners$concordance,3)),cex=1.5)
text(partners$concordance_lb[18]-0.03,18,as.character(round(partners$concordance[18],3)),cex=1.5)
dev.off()

######################################################
# Plot 2: Partners validation set
explorys_all <- read.csv('ccs1_explorys.csv') # Input of plot B

# Sort data by concordance (low to high)
explorys_ehr <- explorys_all[explorys_all$Score_name=='EHR_AF_Final',]
explorys_ehr <- explorys_ehr[order(explorys_ehr$CI_value),]

png('concordance_ccs1_explorys.png',height=500,width=800)
par(oma=c(1,1,1,1))
par(mar=c(3,25,1,1))
col <- ifelse(c((round(mean(explorys_ehr$CI_value),3) >= explorys_ehr$Low_CI) & (round(mean(explorys_ehr$CI_value),3) <= explorys_ehr$High_CI)),'black',
              ifelse(explorys_ehr$Low_CI > mean(explorys_ehr$CI_value),'#6baed6','#e34a33'))

plot(x=explorys_ehr$CI_value,y=seq(1,18,1),xlim=c(0.725,0.875),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col=col,cex=1.5)

axis(1,at=seq(0.675,0.875,0.025),cex.axis=1.6)
axis(2,at=1:18,labels=FALSE,cex.axis=1.6)
mtext('CCS class (level 1)',side=2,line=24.5,cex=1.5)
mtext('Concordance',side=1,line=2.5,cex=1.5)

segments(explorys_ehr$Low_CI,1:18,explorys_ehr$High_CI,1:18,col=col,lwd=2.2)

legend(0.8176,4.19,legend=c('Superior','Expected','Inferior'),lty=1,lwd=3,
       col=c('#6baed6','black','#e34a33'),bty='n',cex=1.5)

plot_names <- explorys_ehr$CCS_PARENT_CATEGORY
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i],' (N=',explorys_ehr$n[i],'; AF=',explorys_ehr$Num_patients_with_AFib[i],')'),
        side=2,line=1,las=2,at=i,cex=1.2)
}

text(explorys_ehr$High_CI[1:17]+0.015,1:17,as.character(round(explorys_ehr$CI_value,3))[1:17],cex=1.5)
text(explorys_ehr$Low_CI[18]-0.015,18,as.character(round(explorys_ehr$CI_value,3))[18],cex=1.5)
dev.off()
