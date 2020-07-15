# Plot c-stat (corrected/uncorrected) by subgroups of age

# Dependencies
library(data.table)

# Load data
compare <- fread(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/cstat_compare.csv')

# Plot
pdf(file='/data/arrhythmia/skhurshid/broad_ibm_afrisk/cstat_compare.pdf',height=3,width=4,pointsize=4)
x <- 1:nrow(compare)
y1 <- compare$c_stat_lb
y2 <- compare$rederived_cstat

plot(x,y1,xaxt='n',xlab='',yaxt='n',ylab='',col='#f46d43',pch=19,bty='n',ylim=c(0.480,0.800))
segments(x,compare$c_stat,x,compare$c_stat_ub,col='#f46d43')

par(new=TRUE)
plot(x,y2,xaxt='n',xlab='',yaxt='n',ylab='',col='#1a9850',pch=19,bty='n',ylim=c(0.480,0.800))
segments(x,compare$rederived_cstat_lb,x,compare$rederived_cstat_ub,col='#1a9850')

axis(1,at=x,labels=compare$age_class)
axis(2,at=seq(0.480,0.800,0.02),las=2)

mtext("age subgroup",1,line=3)
mtext("c-statistic",2,line=3)

legend(x=length(x)-3,y=0.80,c("original c-stat","rederived c-stat"),col=c('#f46d43','#1a9850'),bty='n',pch=19)

dev.off()