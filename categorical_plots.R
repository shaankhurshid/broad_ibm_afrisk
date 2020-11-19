# Plot discrimination in categorical subgroups
pdf(file='/data/arrhythmia/skhurshid/heterogeneity/charge_discrim_cat.pdf',height=4,width=8,pointsize=5)
par(mar=c(6,1,1,1),oma=c(1,1,1,1))

y_variable <- as.numeric(c(as.character(output_sex$c_stat[1:2]),
                           as.character(output_race$c_stat[1:2]),
                           as.character(output_hf$c_stat[1:2]),
                           as.character(output_stroke$c_stat[1:2]),as.character(output_sex$c_stat[3])))
ub <- as.numeric(c(as.character(output_sex$c_stat_ub[1:2]),
                   as.character(output_race$c_stat_ub[1:2]),
                   as.character(output_hf$c_stat_ub[1:2]),
                   as.character(output_stroke$c_stat_ub[1:2]),as.character(output_sex$c_stat_ub[3])))
lb <- as.numeric(c(as.character(output_sex$c_stat_lb[1:2]),
                   as.character(output_race$c_stat_lb[1:2]),
                   as.character(output_hf$c_stat_lb[1:2]),
                   as.character(output_stroke$c_stat_lb[1:2]),as.character(output_sex$c_stat_lb[3])))
n_event <- as.numeric(c(as.character(output_sex$n_event[1:2]),
                        as.character(output_race$n_event[1:2]),
                        as.character(output_hf$n_event[1:2]),
                        as.character(output_stroke$n_event[1:2]),as.character(output_sex$n_event[3])))
col <- c(rep(c('#43a2ca','#e34a33'),4),'darkgray')

plot(x=1:length(y_variable),y=y_variable,col=col,
     bty='n',xaxt='n',yaxt='n',pch=19,xlim=c(0,length(y_variable)),
     ylim=c(0.55,0.85),xlab='',ylab='',cex=2)

axis(1,at=1:length(y_variable),cex.axis=2,
     labels=c('Female','Male','White','Non-white','No HF','HF','No stroke','Stroke','All'))
axis(2,cex.axis=2,at=seq(0.55,0.85,0.05),las=2,pos=0.6)

mtext("Concordance index",2,line=-1,cex=2.5)
mtext("Subgroup",1,line=4,cex=2.5)

text(1:length(y_variable),ub+0.01,as.character(n_event),cex=1.2)

segments(1:length(y_variable),lb,
         1:length(y_variable),ub,
         col=col,lwd=1.5)

dev.off()

# Plot calibration slope in categorical subgroups
pdf(file='/data/arrhythmia/skhurshid/heterogeneity/charge_calib_slope_cat.pdf',height=4,width=8,pointsize=5)
par(mar=c(6,1,1,1),oma=c(1,1,1,1))

y_variable <- as.numeric(c(as.character(output_sex$cal[1:2]),
                           as.character(output_race$cal[1:2]),
                           as.character(output_hf$cal[1:2]),
                           as.character(output_stroke$cal[1:2]),as.character(output_sex$cal[3])))
ub <- as.numeric(c(as.character(output_sex$cal_ub[1:2]),
                   as.character(output_race$cal_ub[1:2]),
                   as.character(output_hf$cal_ub[1:2]),
                   as.character(output_stroke$cal_ub[1:2]),as.character(output_sex$cal_ub[3])))
lb <- as.numeric(c(as.character(output_sex$cal_lb[1:2]),
                   as.character(output_race$cal_lb[1:2]),
                   as.character(output_hf$cal_lb[1:2]),
                   as.character(output_stroke$cal_lb[1:2]),as.character(output_sex$cal_lb[3])))
n_event <- as.numeric(c(as.character(output_sex$n_event[1:2]),
                   as.character(output_race$n_event[1:2]),
                   as.character(output_hf$n_event[1:2]),
                   as.character(output_stroke$n_event[1:2]),as.character(output_sex$n_event[3])))

col <- c(rep(c('#43a2ca','#e34a33'),4),'darkgray')

plot(x=1:length(y_variable),y=y_variable,col=col,
     bty='n',xaxt='n',yaxt='n',pch=19,xlim=c(0,length(y_variable)),
     ylim=c(0,2),xlab='',ylab='',cex=2)

axis(1,at=1:length(y_variable),cex.axis=2,
     labels=c('Female','Male','White','Non-white','No HF','HF','No stroke','Stroke','All'))
axis(2,cex.axis=2,at=seq(0,2,0.25),las=2,pos=0.6)

mtext("Calibration slope",2,line=-1,cex=2.5)
mtext("Subgroup",1,line=4,cex=2.5)

text(1:length(y_variable),ub+0.1,as.character(n_event),cex=1.2)

segments(1:length(y_variable),lb,
         1:length(y_variable),ub,
         col=col,lwd=1.5)

segments(0.6,1,length(y_variable),1,lty=5,col='black')

dev.off()

# Plot relative error in categorical subgroups
pdf(file='/data/arrhythmia/skhurshid/heterogeneity/charge_rel_error_cat.pdf',height=4,width=8,pointsize=5)
par(mar=c(6,1,1,1),oma=c(1,1,1,1))

y_variable <- as.numeric(c(as.character(output_sex$rel_error_mean[1:2]),
                           as.character(output_race$rel_error_mean[1:2]),
                           as.character(output_hf$rel_error_mean[1:2]),
                           as.character(output_stroke$rel_error_mean[1:2]),as.character(output_sex$rel_error_mean[3])))
n_event <- as.numeric(c(as.character(output_sex$n_event[1:2]),
                        as.character(output_race$n_event[1:2]),
                        as.character(output_hf$n_event[1:2]),
                        as.character(output_stroke$n_event[1:2]),as.character(output_sex$n_event[3])))

col <- c(rep(c('#43a2ca','#e34a33'),4),'darkgray')

plot(x=1:length(y_variable),y=y_variable*100,col=col,
     bty='n',xaxt='n',yaxt='n',pch=19,xlim=c(0,length(y_variable)),
     ylim=c(0,60),xlab='',ylab='',cex=2)

axis(1,at=1:length(y_variable),cex.axis=2,
     labels=c('Female','Male','White','Non-white','No HF','HF','No stroke','Stroke','All'))
axis(2,cex.axis=2,at=seq(0,60,10),las=2,pos=0.6)

mtext("Relative prediction error (%)",2,line=-1,cex=2.5)
mtext("Subgroup",1,line=4,cex=2.5)

text(1:length(y_variable),(y_variable*100)+5,as.character(n_event),cex=1.2)

dev.off()