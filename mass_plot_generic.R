# Script for creating density plots for AF risk score distribution in stroke discrim

# Dependencies
library(ggplot2)
library(reshape2)
library(plyr)

# Generate CHADS distribution (controls = no AF, cases = AF)
controls <- count(mydata[is_AFib_combined==0]$pred_risk_CHA2DS2_VASc_Final)
controls_probs <- controls[,2] / nrow(controls)

cases <- count(mydata[is_AFib_combined==1]$pred_risk_CHA2DS2_VASc_Final)
cases_probs <- cases[,2] / nrow(cases)

x <- factor(0:length(cases_probs)) 
y1 <- controls
y2 <- cases

data <- data.frame(x=x,y1=y1,y2=y2)
melted <- melt(data,id="x")

# Plot
png(file='mass_chads.png',height=540,width=760)

ggplot(melted,aes(x=x,y=value,fill=variable)) + geom_bar(stat='identity',position='identity',alpha=0.55,width=1,color='black') +
  scale_x_discrete(breaks=seq(0,9,1),expand=c(0,0),name=expression(paste('   ',CHA[2],DS[2],-VASc,' ',Risk))) +
  scale_y_continuous(breaks=seq(0,0.3,0.05),expand=c(0,0),limits=c(0,0.30),name='frequency') +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('No AF','AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.85,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20))
dev.off()

# Generate C2HEST distribution
controls <- count(mydata[is_AFib_combined==0]$pred_risk_C2HEST_Final)
controls_probs <- controls[,2] / nrow(controls)

cases <- count(mydata[is_AFib_combined==1]$pred_risk_C2HEST_Final)
cases_probs <- cases[,2] / nrow(cases)

x <- factor(0:length(cases_probs)) 
y1 <- controls
y2 <- cases

data <- data.frame(x=x,y1=y1,y2=y2)
melted <- melt(data,id="x")

# Plot
png(file='mass_chest.png',height=540,width=760)

ggplot(melted,aes(x=x,y=value,fill=variable)) + geom_bar(stat='identity',position='identity',alpha=0.55,width=1,color='black') +
  scale_x_discrete(breaks=seq(0,8,1),expand=c(0,0),name=expression(paste('   ',C[2],HEST,' ',Risk))) +
  scale_y_continuous(breaks=seq(0,0.3,0.05),expand=c(0,0),limits=c(0,0.30),name='frequency') +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('No AF','AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.85,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20))
dev.off()