# Dependenices
library(data.table)

# Load output file with hand-curated disease names (OPTIONAL)
hf <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/hf_subgroup.csv')
stroke <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/stroke_subgroup.csv')

# Sort data by concordance (low to high)
output <- output[order(output$concordance),]

png('/data/arrhythmia/skhurshid/broad_ibm_afrisk/concordance_ccs1.png',height=500,width=800)
par(oma=c(1,1,1,1))
par(mar=c(3,32,1,1))
plot(x=output$concordance,y=seq(1,18,1),xlim=c(0.700,1.00),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col='#6baed6',cex=1)
axis(1,at=seq(0.700,1.00,0.025),cex=1.4)
axis(2,at=1:18,labels=FALSE,cex=1.4)
mtext('CCS class (level 1)',side=2,line=31.5,cex=1.3)
mtext('Concordance',side=1,line=2.5,cex=1.3)
segments(output$concordance_lb,1:18,output$concordance_ub,1:18,col='#6baed6',lwd=2.2)

plot_names <- output$ccs1
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i],' (N=',output$n_disease[i],'; AF=',output$incd_af[i],')'),
        side=2,line=1,las=2,at=i)
}

text(output$concordance_ub+0.03,1:17,as.character(round(output$concordance,3)))
text(output$concordance_lb[18]-0.03,18,as.character(round(output$concordance[18],3)))
dev.off()
