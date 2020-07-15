# Dependenices
library(data.table)

# Load subgroup files
f <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/female.csv')
m <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/male.csv')
hf <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/hf.csv')
stroke <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/stroke.csv')
race <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/race.csv')

# Create all
all <- data.frame(Score_Name='EHR_AF_Final',Low_CI=0.782,CI_value=0.783,High_CI=0.784,N_outcome=153151,
                  N_selected_population=4508180,Outcome_inc=NA,N_total=4508180)

# Collapse EHR_AF results only
ehr <- rbind(all[1,],hf[Score_Name=='EHR_AF_Final'],
             stroke[Score_Name=='EHR_AF_Final'],
             race[subgroup=='white',1:8],race[subgroup=='non-white',1:8],
             m[Score_Name=='EHR_AF_Final'],f[Score_Name=='EHR_AF_Final'])
ehr$subgroup <- c('Overall','Heart failure','Stroke','Nonwhite','White','Males','Females')

# Plot
pdf('/data/arrhythmia/skhurshid/broad_ibm_afrisk/subgroup.pdf',height=2,width=4.4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,7.5,1,1))
col <- c('#e31a1c','#ff7f00','#a6d96a','#02818a','#3288bd','#5e4fa2','black')
plot(x=ehr$CI_value,y=seq(1,7,1),xlim=c(0.6,0.825),ylim=c(0.75,7),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=c(18,rep(19,6)),col=rev(col),cex=c(3,rep(1.4,6)),bty='n')
axis(1,at=seq(0.600,0.825,0.025),cex=1.4,pos=0)
axis(2,at=1:7,labels=FALSE,cex=1.4)
mtext('Subgroup',side=2,line=7,cex=1.3)
mtext('Concordance',side=1,line=3.5,cex=1.3)
segments(ehr$Low_CI,1:7,ehr$High_CI,1:7,col=rev(col),lwd=2.2)

plot_names <- c('Overall','Heart failure','Stroke','Nonwhite','White','Male','Female')
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i]),
        side=2,line=1,las=2,at=i)
}

dev.off()
