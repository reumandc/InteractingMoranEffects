#Does the final steps in cleaning. Starting from the results of the more involved cleaning,
#just normalizes kelp time series to an estimate of density per unit usable habitat, and then
#detrends all variables without variance standaridizing.

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): 

#***
#Locations for storing results, and other prep
#***

resloc1<-"../Results/Kelp_DataAfterAllCleaning/"
if (!dir.exists(resloc1))
{
  dir.create(resloc1,recursive=TRUE)
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
#Do the cleaning
#***

#detrend but don't variance standardize the environmental variables
NO3<-wsyn::cleandat(NO3,times,2)$cdat
waves<-wsyn::cleandat(waves,times,2)$cdat

#detrend climate indices
climinds$NPGO<-wsyn::cleandat(climinds$NPGO,times,2)$cdat
climinds$MEI<-wsyn::cleandat(climinds$MEI,times,2)$cdat
climinds$PDO<-wsyn::cleandat(climinds$PDO,times,2)$cdat

#Compute a kelp density per unit useable habitat. To do this, take the max kelp area value ever
#achieved for each habitat patch, use that as the measure of available habitat, and divide the kelp biomass
#time series at that location through by that value. Then detrend.
for (counter in 1:(dim(kelp)[1]))
{
  kelp[counter,]<-kelp[counter,]/(unname(quantile(kelpA[counter,],prob=.9))) #use the 90th percentile for robustness, instead of the max
}
kelpNDT<-kelp
kelp<-wsyn::cleandat(kelp,times,2)$cdat

#***
#Save the final version of cleaned data, save the non-detrended data separately
#***

saveRDS(kelpNDT,file=paste0(resloc1,"Kelp_Quarterly_CleanedFinal_NDT.Rds"))
write.table(kelpNDT,file=paste0(resloc1,"Kelp_Quarterly_CleanedFinal_NDT.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(kelp,file=paste0(resloc1,"Kelp_Quarterly_CleanedFinal.Rds"))
write.table(kelp,file=paste0(resloc1,"Kelp_Quarterly_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(kelpA,file=paste0(resloc1,"KelpA_Quarterly_CleanedFinal.Rds"))
write.table(kelpA,file=paste0(resloc1,"KelpA_Quarterly_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(NO3,file=paste0(resloc1,"NO3_Quarterly_CleanedFinal.Rds"))
write.table(NO3,file=paste0(resloc1,"NO3_Quarterly_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(waves,file=paste0(resloc1,"Waves_Quarterly_CleanedFinal.Rds"))
write.table(waves,file=paste0(resloc1,"Waves_Quarterly_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(locs,file=paste0(resloc1,"Locs_CleanedFinal.Rds"))
write.table(locs,file=paste0(resloc1,"Locs_CleanedFinal.csv"),row.names = FALSE,col.names = TRUE,sep=",")

saveRDS(quarters,file=paste0(resloc1,"Quarters_CleanedFinal.Rds"))
write.table(quarters,file=paste0(resloc1,"Quarters_CleanedFinal.csv"),row.names = FALSE,col.names = TRUE,sep=",")

saveRDS(climinds,file=paste0(resloc1,"Climinds_Quarterly_CleanedFinal.Rds"))
write.table(climinds,file=paste0(resloc1,"Climinds_Quarterly_CleanedFinal.csv"),row.names = FALSE,col.names = TRUE,sep=",")
