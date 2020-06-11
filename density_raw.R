# Script for generating generic density plots for both continuous and discrete scores

# Dependencies
library(ggplot2)
library(rehape2)

setDT(vs)

# Create separate data frames for cases/controls
## Replace *data* with dataset of interest, and replace *af* with disease variable of interest
incident_af <- vs[af_5y_sal==1]
no_incident_af <- vs[af_5y_sal==0]

##### PLOT 3: EHR-AF
# Generate score distribution
x <- list(v1=no_incident_af$score,v2=incident_af$score)
data <- melt(x)

# Density of predicted risk distribution
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) +
  scale_x_continuous(breaks=seq(3,12,1),expand=c(0,0.1),limits=c(3,12)) +   # modify x axis limits as needed
  scale_y_continuous(breaks=seq(0,0.55,0.05),expand=c(0,0),limits=c(0,0.55)) +  # modify y axis limits as needed
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('No AF','AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=18,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=18,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=18),legend.text=element_text(size=18)) +
  labs(x='EHR-AF Score',y='Density')  # modify x axis label as needed (y axis label generally 'density')
ggsave(filename='/Users/Rebecca/Documents/MGH Research/explorys/ehr_density.pdf',height=3,width=3.5,units='in',scale=3,device='pdf')

