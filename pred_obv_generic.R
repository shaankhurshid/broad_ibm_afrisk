# Dependencies
library(survival)

# Script for creating pred/obv plot
## REQUIRES PREDICTED RISK in % AS INPUT

# Create parsing function (this will essentially round predicted risk to the nearest whole number between 0-20)
parse <- function(data,risk){
risk_int <- rep(NA,nrow(data))
data <- data[,risk]
  for (i in 1:20){
    iter <- i-0.5
    for (j in 1:length(data)){
      if(is.na(risk_int[j]) && data[j] < iter){risk_int[j] <- trunc(iter,1)}
    }
    print(paste0('I just finished loop ',i,' out of 20!'))
  }
risk_int <- ifelse(is.na(risk_int),20,risk_int)
return(risk_int)
}

# Apply to scores (REPLACE "RISK" and "X.RISK" with PREDICTED RISKS)
vs$ehraf_pred5_int <- parse(data=vs,risk='risk')
vs$charge_pred5_int <- parse(data=vs,risk='charge.risk')                                                                           
vs$chads_pred5_int <- parse(data=vs,risk='chads.risk')                                                                           
vs$chest_pred5_int <- parse(data=vs,risk='chest.risk')                                                                           

# Function to obtain incidence estimates when given continuous cutoffs (USED FOR EHR-AF AND CHARGE-AF)
observed <- function(data,parsed_risk,time,status){
km <- lower <- upper <- out <- list()
  for (i in 0:(length(unique(data[,parsed_risk]))-1)){
    index <- i+1
    subset <- data[data[,parsed_risk]==i,]
    km[[index]] <- survfit(Surv(subset[,time],subset[,status]) ~ 1, data=subset)$surv
    lower[[index]] <- survfit(Surv(subset[,time],subset[,status]) ~ 1, data=subset)$upper
    upper[[index]] <- survfit(Surv(subset[,time],subset[,status]) ~ 1, data=subset)$lower
    out[[index]] <- c((1-km[[index]][length(km[[index]])])*100,(1-lower[[index]][length(lower[[index]])])*100,(1-upper[[index]][length(upper[[index]])])*100)
    print(paste0('I just finished loop ',i,' out of ',length(unique(data[,parsed_risk]))-1,'!'))
  }
return(out)
}

# Function to obtain incidence estimates when given discrete cutoffs (USED FOR CHADS AND C2HEST)
observed_discrete <- function(data,parsed_risk,time,status){
  km <- lower <- upper <- out <- list()
  index <- 1
  for (i in unique(data[,parsed_risk])[order(unique(data[,parsed_risk]))]){
    subset <- data[data[,parsed_risk]==i,]
    km[[index]] <- survfit(Surv(subset[,time],subset[,status]) ~ 1, data=subset)$surv
    lower[[index]] <- survfit(Surv(subset[,time],subset[,status]) ~ 1, data=subset)$upper
    upper[[index]] <- survfit(Surv(subset[,time],subset[,status]) ~ 1, data=subset)$lower
    out[[index]] <- c((1-km[[index]][length(km[[index]])])*100,(1-lower[[index]][length(lower[[index]])])*100,(1-upper[[index]][length(upper[[index]])])*100)
    print(paste0('I just finished loop ',index,' out of ',length(unique(data[,parsed_risk])),'!'))
    index <- index + 1
  }
  return(out)
}

# Apply to scores (DATA = DATA, PARSED_RISK = OUTPUT OF PARSE FUNCTION, TIME=EVENT TIME, STATUS=AF YES/NO)
ehraf_obv <- do.call('rbind',observed(data=vs,parsed_risk='ehraf_pred5_int',time='af_5y_sal.t',status='af_5y_sal'))
charge_obv <- do.call('rbind',observed(data=vs,parsed_risk='charge_pred5_int',time='af_5y_sal.t',status='af_5y_sal'))
chads_obv <- do.call('rbind',observed_discrete(data=vs,parsed_risk='chads_pred5_int',time='af_5y_sal.t',status='af_5y_sal'))
chest_obv <- do.call('rbind',observed_discrete(data=vs,parsed_risk='chest_pred5_int',time='af_5y_sal.t',status='af_5y_sal'))

# Plotting
# Open file connection
png("/data/arrhythmia/skhurshid/ehr_af/pred_obv2.png",height=760,width=1200)

## Params
par(oma=c(1,1,0,0))

## x-scales
x2 <- -0.375:19.625
x3 <- -0.125:19.875
x1 <- unique(vs$chads_pred5_int)[order(unique(vs$chads_pred5_int))]+0.125
x4 <- unique(vs$chest_pred5_int)[order(unique(vs$chest_pred5_int))]+0.375

## colors
col1 <- "#1f78b4"; col2 <- "#33a02c"; col3 <- "#e31a1c"; col4 <- "#ff7f00"

# Plot
## Point estimates (YOU MAY NEED TO TWEAK THIS BASED ON HOW MANY POINTS CHADS AND C2HEST HAVE [THIS IS UNIQUE TO DATASET])
plot(x=c(x1,x2,x3,x4),y=c(chads_obv[,1],charge_obv[,1],ehraf_obv[,1],chest_obv[,1]),
     pch=19,col=c(rep(col1,5),rep(col2,21),rep(col3,21),rep(col4,4)),xlim=c(-0.5,21),ylim=c(0,25),
     xaxt='n',yaxt='n',xlab='',ylab='',cex=1.4,frame.plot = F)

## CIs
segments(x3,ehraf_obv[,2],x3,ehraf_obv[,3],col=col3)
segments(x2,charge_obv[,2],x2,charge_obv[,3],col=col2)
segments(x1,chads_obv[,2],x1,chads_obv[,3],col=col1)
segments(x4,chest_obv[,2],x4,chest_obv[,3],col=col4)

## Axes
axis(side=1,cex.axis=1.5,at=seq(0,20,1),labels=seq(0,20,1),pos=-0.5)
axis(side=2,cex.axis=1.5,at=seq(0,25,5),labels=seq(0,25,5),las=1,pos=-0.5)

## Labels
mtext(side=2,"Observed 5-Year AF (%)",line=1.8,cex=1.5)
mtext(side=1,"Predicted 5-Year AF Risk (%)",line=3,cex=1.5)

## Segments
segments(0,-0.5,-0.5,-0.5)
segments(-0.5,0,-0.5,-0.5)
segments(20,-0.5,20.5,-0.5)

# 45 Line
black_trans <- adjustcolor('black', alpha.f = 0.6) 
segments(-0.5,-0.5,20,20,col=black_trans)

## Legend
legend(-0.2,24.5,c('CHADS-VASc','CHARGE-AF','EHR-AF','C2HEST'),
       pch=19,col=c(col1,col2,col3,col4),cex=1.3,bty='n')

dev.off()
