# Dependenices
library(data.table)

# Load subgroup files
f <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/female.csv')
m <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/male.csv')
hf <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/hf.csv')
stroke <- fread('/data/arrhythmia/skhurshid/broad_ibm_afrisk/stroke.csv')

# Create all
all <- data.frame(Score_Name='EHR_AF_Final',Low_CI=0.780,CI_value=0.781,High_CI=0.782,N_outcome=152205,
                  N_selected_population=4510814,Outcome_inc=NA,N_total=4510814)

# Collapse EHR_AF results only
ehr <- rbind(f[Score_Name=='EHR_AF_Final'],m[Score_Name=='EHR_AF_Final'],
             stroke[Score_Name=='EHR_AF_Final'],hf[Score_Name=='EHR_AF_Final'],all[1,])
ehr$subgroup <- c('Females','Males','Stroke','Heart failure','Overall')

# Plot
png('/data/arrhythmia/skhurshid/broad_ibm_afrisk/subgroup.png',height=350,width=800,res=100)
par(oma=c(1,1,1,1))
par(mar=c(3,19,1,1))
col <- c('#e31a1c','#1f78b4','#33a02c','#ff7f00','black')
plot(x=ehr$CI_value,y=seq(5,1,-1),xlim=c(0.600,0.820),ylim=c(0.75,5.25),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col=col,cex=1.4)
axis(1,at=seq(0.600,0.820,0.020),cex=1.4)
axis(2,at=1:18,labels=FALSE,cex=1.4)
mtext('Subgroup',side=2,line=18.5,cex=1.3)
mtext('Concordance',side=1,line=2.5,cex=1.3)
segments(ehr$Low_CI,5:1,ehr$High_CI,5:1,col=col,lwd=2.2)

plot_names <- c('overall','heart failure','stroke','male','female')
for (i in 1:length(plot_names)){
  mtext(paste0(plot_names[i],' (N=',ehr$N_selected_population[i],'; AF=',ehr$N_outcome[i],')'),
        side=2,line=1,las=2,at=i)
}

dev.off()
