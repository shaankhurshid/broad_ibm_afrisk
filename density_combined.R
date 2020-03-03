# Script for generating generic density plots for both continuous and discrete scores

# Dependencies
library(ggplot2)
library(rehape2)

# Create separate data frames for cases/controls
## Replace *data* with dataset of interest, and replace *af* with disease variable of interest
incident_af <- mydata[is_AFib_combined==1,]
no_incident_af <- mydata[is_AFib_combined==0,]

##### PLOT 1: CHADS
controls <- count(no_incident_af$pred_risk_CHA2DS2_VASc_Final)
controls_probs <- controls[,2] / nrow(no_incident_af)

cases <- count(incident_af$pred_risk_CHA2DS2_VASc_Final)
cases_probs <- cases[,2] / nrow(incident_af)

x <- factor(cases[,1])
y1 <- controls_probs
y2 <- cases_probs

data <- data.frame(x=x,y1=y1,y2=y2)
melted <- melt(data,id="x")

# Plot
ggplot(melted,aes(x=x,y=value,fill=variable)) + geom_bar(stat='identity',position='identity',alpha=0.55,width=1,color='black') +
  scale_x_discrete(breaks=x,expand=c(0,0),name=expression(paste('   ',CHA[2],DS[2],-VASc,' ',Risk))) +
  scale_y_continuous(breaks=seq(0,0.3,0.05),expand=c(0,0),limits=c(0,0.30),name='frequency') +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('No AF','AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.85,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20))
ggsave(filename='chads_density.pdf',height=3,width=3.5,units='in',scale=4,device='pdf')

##### PLOT 2: CHEST
controls <- count(no_incident_af$pred_risk_C2HEST_Final)
controls_probs <- controls[,2] / nrow(controls)

cases <- count(incident_af$pred_risk_C2HEST_Final)
cases_probs <- cases[,2] / nrow(cases)

x <- factor(cases[,1])
y1 <- controls_probs
y2 <- cases_probs

data <- data.frame(x=x,y1=y1,y2=y2)
melted <- melt(data,id="x")

# Plot
ggplot(melted,aes(x=x,y=value,fill=variable)) + geom_bar(stat='identity',position='identity',alpha=0.55,width=1,color='black') +
  scale_x_discrete(breaks=x,expand=c(0,0),name=expression(paste('   ',C[2],HEST,' ',Risk))) +
  scale_y_continuous(breaks=seq(0,0.3,0.05),expand=c(0,0),limits=c(0,0.30),name='frequency') +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('No AF','AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.85,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20))
ggsave(filename='chest_density.pdf',height=3,width=3.5,units='in',scale=4,device='pdf')

##### PLOT 3: CHARGE-AF
# Generate score distribution
x <- list(v1=incident_af$pred_risk_CHARGE_Final,v2=no_incident_af$pred_risk_CHARGE_Final)
data <- melt(x)

# Density of predicted risk distribution
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) +
  scale_x_continuous(breaks=seq(0,15,1),expand=c(0,0.1),limits=c(0,15)) +   # modify x axis limits as needed
  scale_y_continuous(breaks=seq(0,0.90,0.05),expand=c(0,0),limits=c(0,0.90)) +  # modify y axis limits as needed
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('incident AF','no incident AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted 5-year AF risk using CHARGE-AF score (%)',y='density')  # modify x axis label as needed (y axis label generally 'density')
ggsave(filename='charge_density.pdf',height=3,width=3.5,units='in',scale=4,device='pdf')

##### PLOT 3: EHR-AF
# Generate score distribution
x <- list(v1=incident_af$pred_risk_EHR_Final,v2=no_incident_af$pred_risk_EHR_Final)
data <- melt(x)

# Density of predicted risk distribution
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) +
  scale_x_continuous(breaks=seq(0,15,1),expand=c(0,0.1),limits=c(0,15)) +   # modify x axis limits as needed
  scale_y_continuous(breaks=seq(0,0.90,0.05),expand=c(0,0),limits=c(0,0.90)) +  # modify y axis limits as needed
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('incident AF','no incident AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted 5-year AF risk using EHR-AF score (%)',y='density')  # modify x axis label as needed (y axis label generally 'density')
ggsave(filename='ehr_density.pdf',height=3,width=3.5,units='in',scale=4,device='pdf')

