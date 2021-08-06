#Just plots time series after the "more involved" cleaning (but prior to the final cleaning done in 
#AnalysisStep_LinearModels.R), for visual inspection.

rm(list=ls())

#***
#Locations for storing results
#***

resloc<-"../Results/Kelp_TimeSeriesPlots/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***
#Load the data
#***

#load the data
datloc<-"../Results/Kelp_DataAfterMoreInvolvedCleaning/"
kelp<-readRDS(paste0(datloc,"Kelp_Quarterly_CleanedMoreInvolved.Rds")) 
kelpA<-readRDS(paste0(datloc,"KelpA_Quarterly_CleanedMoreInvolved.Rds")) 
NO3<-readRDS(paste0(datloc,"NO3_Quarterly_CleanedMoreInvolved.Rds")) 
waves<-readRDS(paste0(datloc,"Waves_Quarterly_CleanedMoreInvolved.Rds")) 
locs<-readRDS(paste0(datloc,"Locs_CleanedMoreInvolved.Rds"))
quarters<-readRDS(paste0(datloc,"Quarters_CleanedMoreInvolved.Rds"))
climinds<-readRDS(paste0(datloc,"Climinds_Quarterly_CleanedMoreInvolved.Rds"))
times<-1:(dim(quarters)[1])

numts<-dim(kelp)[1]
lents<-dim(kelp)[2]

#***
#Make the plots
#***

pan.wd<-3
pan.ht<-pan.wd
gap<-0.2
xax.ht<-.6
yax.wd<-xax.ht
tot.wd<-5*yax.wd+5*gap+5*pan.wd
tot.ht<-xax.ht+pan.ht+gap

quat<-0:(dim(kelp)[2]-1)

for (counter in 1:(dim(kelp)[1]))
{
  #plot the quarterly data
  png(paste0(resloc,"SitePlot_",counter,".png"),res=300,units="in",width = tot.wd,height = tot.ht)
  
  #plot the kelp biomass on the first set of axes
  par(fig=c(yax.wd/tot.wd,
            (yax.wd+pan.wd)/tot.wd,
            (xax.ht)/tot.ht,
            (xax.ht+pan.ht)/tot.ht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0))
  plot(quat,kelp[counter,],type="l",col="green",xlab="quarter",ylab="kelp") 
  mtext("quarter",side=1,line=1.6)
  mtext("kelp",side=2,line=1.6)
  
  #kelp area on the second axes
  par(fig=c((2*yax.wd+pan.wd+gap)/tot.wd,
            (2*yax.wd+2*pan.wd+gap)/tot.wd,
            (xax.ht)/tot.ht,
            (xax.ht+pan.ht)/tot.ht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  plot(quat,kelpA[counter,],type="l",col="green",xlab="quarter",ylab="kelp") 
  mtext("quarter",side=1,line=1.6)
  mtext("kelp area",side=2,line=1.6)
  
  #plot nitrates on the next set of axes
  par(fig=c((3*yax.wd+2*pan.wd+2*gap)/tot.wd,
            (3*yax.wd+3*pan.wd+2*gap)/tot.wd,
            (xax.ht)/tot.ht,
            (xax.ht+pan.ht)/tot.ht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  plot(quat,NO3[counter,],type="l",col="black",xlab="quarter",ylab="NO3") 
  mtext("quarter",side=1,line=1.6)
  mtext("NO3",side=2,line=1.6)
  
  #plot waves on the next set of axes
  par(fig=c((4*yax.wd+3*pan.wd+3*gap)/tot.wd,
            (4*yax.wd+4*pan.wd+3*gap)/tot.wd,
            (xax.ht)/tot.ht,
            (xax.ht+pan.ht)/tot.ht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  plot(quat,waves[counter,],type="l",col="blue",xlab="quarter",ylab="waves") 
  mtext("quarter",side=1,line=1.6)
  mtext("waves",side=2,line=1.6)
  
  #plot all the sites, with this one identified
  par(fig=c((5*yax.wd+4*pan.wd+4*gap)/tot.wd,
            (5*yax.wd+5*pan.wd+4*gap)/tot.wd,
            (xax.ht)/tot.ht,
            (xax.ht+pan.ht)/tot.ht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  plot(locs$Lon,locs$Lat,type="p",pch=20,cex=.1,xlab="Lon",ylab="Lat")
  points(locs$Lon[counter],locs$Lat[counter],type="p",col="red",pch=20,cex=.5)
  mtext("Lon",side=1,line=1.6)
  mtext("Lat",side=2,line=1.6)
  
  dev.off()
}
