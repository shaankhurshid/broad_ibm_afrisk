# Script to visualize calibration by decile

# Quantile sorter
classifier <- function(risk,ncuts){
  cuts <- quantile(risk,probs=seq(0,1,1/ncuts))
  index <- rep(NA,length(risk))
  for (i in 1:(length(cuts)-1)){
    for (j in 1:length(risk)){
      index[j] <- ifelse(risk[j] >= cuts[i],i,index[j])}}
  return(index)
}

# Use classifier to classify scores into quantiles (size 10)
fhs_hf_m$fhs_hf_quantile <- classifier(risk=fhs_hf_m$hf_pred4,ncuts=20)
fhs_hf_f$fhs_hf_quantile <- classifier(risk=fhs_hf_f$hf_pred4,ncuts=20)
fhs_hf_m$fhs_hf_quantile_cal <- classifier(risk=fhs_hf_m$hf_pred4_cal,ncuts=20)
fhs_hf_f$fhs_hf_quantile_cal <- classifier(risk=fhs_hf_f$hf_pred4_cal,ncuts=20)

### CALCULATE OBSERVED RISK IN EACH QUANTILE
# Observed risks
fhs_hf_obv_m <- fhs_hf_m[,lapply(.SD,mean),by='fhs_hf_quantile',.SDcols='inpt_hf_4y'][order(fhs_hf_quantile),'inpt_hf_4y']
fhs_hf_obv_f <- fhs_hf_f[,lapply(.SD,mean),by='fhs_hf_quantile',.SDcols='inpt_hf_4y'][order(fhs_hf_quantile),'inpt_hf_4y']
fhs_hf_obv_m_cal <- fhs_hf_m[,lapply(.SD,mean),by='fhs_hf_quantile_cal',.SDcols='inpt_hf_4y'][order(fhs_hf_quantile_cal),'inpt_hf_4y']
fhs_hf_obv_f_cal <- fhs_hf_f[,lapply(.SD,mean),by='fhs_hf_quantile_cal',.SDcols='inpt_hf_4y'][order(fhs_hf_quantile_cal),'inpt_hf_4y']

### CALCULATE AVERAGE PREDICTED RISK IN EACH QUANTILE
fhs_hf_pred_m <- fhs_hf_m[,mean(hf_pred4),by="fhs_hf_quantile"][order(fhs_hf_quantile)]
fhs_hf_pred_f <- fhs_hf_f[,mean(hf_pred4),by="fhs_hf_quantile"][order(fhs_hf_quantile)]
fhs_hf_pred_m_cal <- fhs_hf_m[,mean(hf_pred4_cal),by="fhs_hf_quantile_cal"][order(fhs_hf_quantile_cal)]
fhs_hf_pred_f_cal <- fhs_hf_f[,mean(hf_pred4_cal),by="fhs_hf_quantile_cal"][order(fhs_hf_quantile_cal)]

# Plots for visualization
## Male
pdf(file='/data/arrhythmia/skhurshid/hf_prediction/fhs_hf_cal_qual_m_recalibrated.pdf',height=4,width=8,
    pointsize=3)
par(oma=c(2,3,1,1))
par(oma=c(2,3,1,1))

x <- fhs_hf_pred_m_cal$V1
y <- do.call(rbind,fhs_hf_obv_m_cal)[1,]*100

plot(x,y,xlab='',ylab='',yaxt='n',
     xaxt='n',xlim=c(0,20),ylim=c(0,20),pch=19,cex=1.5)

axis(1,at=seq(0,20,5),cex.axis=1.7)
axis(2,cex.axis=1.6,at=seq(0,20,5),las=1)

segments(-1,-1,31,31,lwd=1.2,lty=2)

mtext("Men",side=3,cex=3,line=2,at=13.5)
mtext("Predicted risk of HF at 4 years (%)",side=1,cex=1.7,line=3)
mtext("Incidence of HF at 4 years (%)",side=2,cex=1.7,line=4.5)

dev.off()
