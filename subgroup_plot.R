# Dependenices
library(data.table)

# Load subgroup files
f <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/female.csv')
m <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/male.csv')
hf <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/hf.csv')
stroke <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/stroke.csv')

# Create all
all <- data.frame(Score_Name='EHR_AF_Final',Low_CI=0.782,CI_value=0.783,High_CI=0.784,N_outcome=153151,
                  N_selected_population=4508180,Outcome_inc=NA,N_total=4508180)

# Collapse EHR_AF results only
ehr <- rbind(all[1,],hf[Score_Name=='EHR_AF_Final'],
             stroke[Score_Name=='EHR_AF_Final'],
             m[Score_Name=='EHR_AF_Final'],f[Score_Name=='EHR_AF_Final'])
ehr$subgroup <- c('Overall','Heart failure','Stroke','Males','Females')

# Plot
pdf('/data/arrhythmia/skhurshid/broad_ibm_afrisk/subgroup.pdf',height=1.5,width=4.4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,7,1,1))
col <- c('black','#1f78b4','#33a02c','#ff7f00','#e31a1c')
plot(x=ehr$CI_value,y=seq(1,5,1),xlim=c(0.6,0.85),ylim=c(0.75,5),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=c(18,rep(19,4)),col=col,cex=c(3,rep(1.4,4)),bty='n')
axis(1,at=seq(0.6,0.85,0.025),cex=1.4,pos=0)
axis(2,at=1:5,labels=FALSE,cex=1.4)
mtext('Subgroup',side=2,line=6.5,cex=1.3)
mtext('Concordance',side=1,line=3.5,cex=1.3)
segments(ehr$Low_CI,1:5,ehr$High_CI,1:5,col=col,lwd=2.2)

plot_names <- c('Overall','Heart failure','Stroke','Male','Female')
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i]),
        side=2,line=1,las=2,at=i)
}

dev.off()
