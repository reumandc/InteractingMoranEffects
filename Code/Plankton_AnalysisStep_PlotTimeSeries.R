#Just plots time series after the "basic" cleaning (but prior to the final cleaning done in 
#Plankton_AnalysisStep_LinearModels.R), for visual inspection.

rm(list=ls())

#***
#Locations for storing results
#***

resloc<-"../Results/Plankton_TimeSeriesPlots/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***
#Load the data
#***

#load the data
datloc<-"../Results/Plankton_DataAfterBasicCleaning/"
pci<-readRDS(paste0(datloc,"Pci_CleanedBasic.Rds")) 
calfin<-readRDS(paste0(datloc,"Calfin_CleanedBasic.Rds")) 
temp<-readRDS(paste0(datloc,"Temp_CleanedBasic.Rds")) 
year<-readRDS(paste0(datloc,"Year_CleanedBasic.Rds"))
site<-readRDS(paste0(datloc,"Site_CleanedBasic.Rds"))

numts<-dim(pci)[1]
lents<-dim(pci)[2]

#***
#Make the plots
#***

pan.wd<-3
pan.ht<-pan.wd
gap<-0.2
xax.ht<-.6
yax.wd<-xax.ht
tot.wd<-3*yax.wd+3*gap+3*pan.wd
tot.ht<-xax.ht+pan.ht+gap

for (counter in 1:(dim(pci)[1]))
{
  #plot the pci data
  png(paste0(resloc,"SitePlot_",counter,".png"),res=300,units="in",width = tot.wd,height = tot.ht)
  
  #plot pci on the first set of axes
  par(fig=c(yax.wd/tot.wd,
            (yax.wd+pan.wd)/tot.wd,
            (xax.ht)/tot.ht,
            (xax.ht+pan.ht)/tot.ht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0))
  plot(year,pci[counter,],type="l",col="green",xlab="year",ylab="pci") 
  mtext("year",side=1,line=1.6)
  mtext("pci",side=2,line=1.6)
  
  #temp on the second axes
  par(fig=c((2*yax.wd+pan.wd+gap)/tot.wd,
            (2*yax.wd+2*pan.wd+gap)/tot.wd,
            (xax.ht)/tot.ht,
            (xax.ht+pan.ht)/tot.ht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  plot(year,temp[counter,],type="l",col="red",xlab="year",ylab="temperature") 
  mtext("year",side=1,line=1.6)
  mtext("temperature",side=2,line=1.6)
  
  #plot cal fin on the next set of axes
  par(fig=c((3*yax.wd+2*pan.wd+2*gap)/tot.wd,
            (3*yax.wd+3*pan.wd+2*gap)/tot.wd,
            (xax.ht)/tot.ht,
            (xax.ht+pan.ht)/tot.ht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  plot(year,calfin[counter,],type="l",col="blue",xlab="year",ylab="cal fin") 
  mtext("year",side=1,line=1.6)
  mtext("cal fin",side=2,line=1.6)
  
  dev.off()
}

