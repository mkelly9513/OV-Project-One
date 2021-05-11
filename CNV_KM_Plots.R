#the data needs to be formatted to match a certain input 
#in this case the format is as follows 
#col 1: patient id 
#col 2: Surival event, 1 being death
#col 3: Survival Time
#col 4: patient id
#col 5: first feature
#cols 6-end: rest of the features 

#what you need is survival events, survival times, and feature information as the columns 
#all of the patients/samples are the rows 


library(survival)
Surv<- survival::Surv


###########################################################
##########  Function Definition               #############
###########################################################

kmPlot<-function(x,
                 event,
                 stime,
                 varName="",
                 ymin=0,
                 lineColors=NA,
                 nclasses=NA,
                 ci=F,
                 ylab="Overall Survival (Probability)",
                 xlab="Months",
                 plegloc="bottomright",
                 lwid=1,
                 main.text.size=1,
                 citime=c(5,10)){
  
  mainLabel <- varName
  event<-as.numeric(as.vector(event))
  stime<-as.numeric(as.vector(stime))
  if(is.numeric(x)){
    tclass <- factor(x[!is.na(x)])
    event <- event[!is.na(x)]
    stime <- stime[!is.na(x)]
  }else{
    tclass <- factor(x[x!=""])
    event <- event[x!=""]
    stime <- stime[x!=""]
  }
  
  if(is.na(nclasses)){
    nclasses<-nlevels(tclass)
  }
  
  if(length(lineColors)<=1){
    lineColors<-seq(1,nclasses)
  }
  
  y<-survfit(Surv(stime, event)~tclass)
  
  plot(y,
       col=lineColors,
       main=mainLabel,
       ylim=c(ymin,1),
       ylab=ylab,
       xlab=xlab,
       cex.main=main.text.size,
       lwd=lwid)
  
  sy<-summary(y)
  
  if(ci==T)
  {
    for(j in 2:(nclasses+1))
    {
      thislevel<-which(sy$strata==levels(sy$strata)[j-1])
      if(levels(sy$strata)[j-1] != "Basal")
      {
        for(k in 1:length(citime))
        {
          closest<-min(abs(citime[k]-sy$time[thislevel]))==abs(citime[k]-sy$time[thislevel])
          low<-sy$lower[thislevel][closest]
          high<-sy$upper[thislevel][closest]
          segments(citime[k],low,citime[k],high,lty=2)
          segments(citime[k]-0.3,low,citime[k]+0.3,low)
          segments(citime[k]-0.3,high,citime[k]+0.3,high)
          print(paste(levels(sy$strata)[j-1],k,low,sy$surv[thislevel][closest],high))
        }
      }
    }
  }
  legendtxt<- paste(levels(tclass)," ",table(tclass,event)[,2],"/",y$n,sep="")
  legend(0,0.2,legend=legendtxt,col=lineColors,bty="n",lty=rep(1,nclasses),lwd=lwid)
  
  #joelsway
  pvalue<-1-pchisq(survdiff(Surv(stime, event)~tclass)$chisq,nclasses-1)
  
  if(pvalue<0.05)
  {
    legend(plegloc,legend=paste("log rank p = ",signif(pvalue,3),sep=""),cex=1.2,text.col="red",bty="n")
  }else{
    legend(plegloc,legend=paste("log rank p=",signif(pvalue,3),sep=""),cex=1.2,bty="n")
  }	
}

##############################################################################

#setwd("~/Documents/KaplanMPlots/OVCAR_Curated/")
##############################################################################
########################  PDF KM Plots  ########################
##############################################################################
#CNV data 
#replace with your file 
km.data<-read.delim("OVCAR_Tumor_Survival_With_CNV_Data_cleaned_no2.txt",sep="\t",header=T)


#tell it what the survival column is 
survival.time <- km.data[,3]     # Time_5yr
#tell it the event 
survival.event <- km.data[,2]    # Event_5yr
# if you want to do binary things you can use feat 0 and replace the features.1 with features.0 in the code 
features.0 <- 5 
#all SE windows 336 total windows  
features.1 <- 5:340
#se match both 
#features.1 <- 126:128

#oneregion
#features.1 <- 128

#for running of code we clip off the first patient ID
km.data<-km.data[,-1]

#set the name of the PDF/s
pdf("OVCAR_SEtest_Cleaned.KMPLOTS.Analysis.pdf", width=12, height=7)
par(mgp=c(1.3,.35,.0), mai=c(.5,.5,.4,.2), mfcol=c(1,2))

#km.data<- km.data.2
#clean it to remove the NAs for analysis 
km.data<-km.data[is.na(survival.event) == FALSE,]
survival.time<- as.numeric(survival.time[is.na(survival.event) == FALSE])
survival.event<- as.numeric(survival.event[is.na(survival.event) == FALSE])
for (i in features.1)
{
  
  quantile.00<-as.numeric(quantile(km.data[,i], 0.00))
  quantile.33<-as.numeric(quantile(km.data[,i], 0.3333))
  quantile.50<-as.numeric(quantile(km.data[,i], 0.50))
  quantile.66<-as.numeric(quantile(km.data[,i], 0.6666))
  #quantile.100<-as.numeric(quantile(km.data[,i], 0.100))
  quantile.100<-as.numeric(quantile(km.data[,i], 1))
  
  #G.2<- cut(km.data[,i],c(-1000,quantile.50,1000),c("low","high"))
  #this one is for rna seq G.2<- cut(km.data[,i],c(quantile.00,quantile.50,quantile.100),c("low","high"))
  G.2<- cut(km.data[,i],c(quantile.00,quantile.50,quantile.100),c("low","high"))
  #G.2<- cut(km.data[,i],c(-2,0,2),c('del','amp'))
  #G.3<- cut(km.data[,i],c(-1000,quantile.33,quantile.66,1000),c("low","med","high"))
  G.3<- cut(km.data[,i],c(quantile.00,quantile.33,quantile.66,quantile.100),c("low","med","high"))
  
  ################## Group 2 ######################################################
  kmPlot(G.2, 
         survival.event, 
         survival.time, 
         ylab="Event Free Survival (Probability)", #"Overall Survival (Probability)", "Relapse Free Survival (Probability)"
         xlab="Days", # "Years", "Months"
         lwid=1,
         lineColors=c("green","red"), 
         main.text.size=0.8,
         varName=names(km.data)[i])
  
  ################## Group 3 ######################################################
  kmPlot(G.3, 
         survival.event, 
         survival.time, 
         ylab="Event Free Survival (Probability)", #"Overall Survival (Probability)", "Relapse Free Survival (Probability)"
         xlab="Days", # "Years", "Months"
         lwid=1,
         lineColors=c("green","blue","red"), 
         main.text.size=0.8,
         varName=names(km.data)[i])
}
dev.off()

#coxph 

coxph(formula = Surv(survival.time, survival.event) ~ G.2 )


