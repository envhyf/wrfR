#'  @title Basic statistical estimators
#'  @description
#'  \code{basic.stats} is a simple basic statistical estimators calculation
#'  @param data vector
#'  @return data.frame
#'  @author I. Lopez-Coto, 2012 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export

basic.stats<-function(data)
{

  data.mean <- mean(data)
  data.max  <- max(data)
  data.min  <- min(data)
  data.sd   <- sd(data)
  
    data.quantile.95<-quantile(data,0.95,na.rm=TRUE)
    data.quantile.05<-quantile(data,0.05,na.rm=TRUE)
    data.quantile.50<-quantile(data,0.5,na.rm=TRUE)
  
  results<-data.frame(data.mean,data.sd,data.min,data.quantile.05,data.quantile.50,data.quantile.95,data.max)
  dimnames(results)[[2]]<-c('mean','sd','min','P05','P50','P95','max')
  dimnames(results)[[1]]<- ''
  
  return(results)
  
}

#'  @title Verification methods
#'  @description
#'  \code{VM} performes some verification methods
#'  @param data.obs vector of observations
#'  @param data.fcst vector of forecast data
#'  @param plotMode \code{logical}
#'  @param main string to be used on the plots title 
#'  @return data.frame
#'  @author I. Lopez-Coto, 2012 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export

VM <- function(data.obs,data.fcst,plotMode=TRUE,main_info=NULL)
{
require(verification)

delta <- data.fcst-data.obs
bias <- mean(delta)
rmsd <- sqrt(mean((delta)^2))

COR<-cor(data.fcst,data.obs,method = "pearson")
CORk<-cor(data.fcst,data.obs,method = "kendall")

stderr <- sqrt(mean((delta-bias)^2))

LEPS<-leps(data.fcst,data.obs,plot = FALSE)$leps.0  #Linear Error in Probability Space (better close to 0)
lmrs<-lm(data.fcst ~ data.obs)
SLM<-lmrs$coefficients[2]
R2<-summary(lmrs)$r.squared

LMR<-data.frame(bias,rmsd,stderr,SLM,R2,LEPS,COR,CORk)

if (plotMode){
# Setup a 2x2 plotting page
# x11()
old.par<-par()
par(mfrow=c(2,2))
par(bg = "white")
# Create histograms of the forecsat and observation values
chist <- hist(c(data.fcst, data.obs), plot=FALSE)
brks  <- chist$breaks
fhist <- hist(data.fcst, breaks=brks, plot=FALSE)
ohist <- hist(data.obs,  breaks=brks, plot=FALSE)
ymax  <- max(fhist$counts, ohist$counts)


title <- paste("Forecast Histogram for", main_info)
hist(data.fcst, main=title, ylab="Frequency", xlab="Forecast",
     breaks=brks, ylim=c(0, ymax))
title <- paste("Observation Histogram for ", main_info)
hist(data.obs, main=title, ylab="Frequency", xlab="Observation",
     breaks=brks, ylim=c(0, ymax))

# Colors for plotting points
colors <- rgb(0, 0, 1, 0.25)

# Create a scatter plot
title <- paste("Scatter Plot for ", main_info,"\n","Slope = ",formatC(LMR$SLM,digits=3),"; R2 = ",formatC(LMR$R2,digits=3))
plot(y=data.fcst, x=data.obs, main=title, xlab="Observation", ylab="Forecast",
     col=colors, pch=19)
abline(a=0, b=1, lwd=2, lty=2)

# Create a Q-Q plot
title <- paste("Q-Q Plot for ", main_info, "\n","LEPS = ",formatC(LMR$LEPS,digits=3))
qqplot(y=data.fcst, x=data.obs, main=title, xlab="Observation", ylab="Forecast", col=colors, pch=19)
abline(a=0, b=1, lwd=2, lty=2)

#dev.copy(png,paste("LMR_",main_info,".png",sep=""),height=540,width=640)
#dev.off()
#graphics.off()
par(old.par)
}

return(LMR)
}
