# Script for generating generic density plots for two conditions

# Dependencies
library(ggplot2)

# Create separate data frames for cases/controls
## Replace 'data' with dataset of interest, and replace 'af' with disease variable
incident_af <- data[af==1,]
no_incident_af <- data[af==0,]

### Density Plots
# Generate CHARGE distribution
# Replace 'predicted_risk' with predicted risk marker
x <- list(v1=incident_af$predicted_risk,v2=no_incident_af$predicted_risk)
data <- melt(x)

# Density of predicted risk distribution
png(file='density_plot.png',height=540,width=760)
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) +
  scale_x_continuous(breaks=seq(0,15,1),expand=c(0,0.1),limits=c(0,15)) +
  scale_y_continuous(breaks=seq(0,0.90,0.05),expand=c(0,0),limits=c(0,0.90)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('incident AF','no incident AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted 5-year AF risk using XXX score (%)',y='density')
dev.off()
