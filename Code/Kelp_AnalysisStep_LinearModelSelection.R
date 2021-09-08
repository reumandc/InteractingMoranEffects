#Does some linear models and information criterion analyses to select numbers of lags to use in 
#ARMA-type models of kelp.

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): ncf

source("PickLags.R")

#***
#Locations for storing results, and other prep
#***

resloc<-"../Results/Kelp_LinearModelSelection/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***
#Load the data
#***

#load the data
datloc<-"../Results/Kelp_DataAfterAllCleaning/"
kelp<-readRDS(paste0(datloc,"Kelp_Quarterly_CleanedFinal.Rds")) 
NO3<-readRDS(paste0(datloc,"NO3_Quarterly_CleanedFinal.Rds")) 
waves<-readRDS(paste0(datloc,"Waves_Quarterly_CleanedFinal.Rds")) 
locs<-readRDS(paste0(datloc,"Locs_CleanedFinal.Rds"))
quarters<-readRDS(paste0(datloc,"Quarters_CleanedFinal.Rds"))
climinds<-readRDS(paste0(datloc,"Climinds_Quarterly_CleanedFinal.Rds"))
times<-1:(dim(quarters)[1])

numts<-dim(kelp)[1]
lents<-dim(kelp)[2]

#***
#get your regions of CA, copied from another script
#***

#from google maps
MontBLatLon<-c(36.62393626182818, -121.8511283749223) #Monterey Bay, closer to the southern end
CarBLatLon<-c(36.52316357182972, -121.95535513873469) #The point to the south of Carmel Bay
MorBLatLon<-c(35.427553615598676, -120.89133343717805) #Morro Bay, closer to the northern end
PtConcLatLon<-c(34.44830424836278,-120.47125509005721) #Pt Conception
SBLatLon<-c(34.407876790966135, -119.68428533917637) #Santa Barbara pier
OxLatLon<-c(34.2310543912163, -119.2726684255394) #Oxnard, where the Santa Clara River empties
LALatLon<-c(34.00907538429378, -118.50422514984244) #LA, near the Santa Monica pier

dists<-ncf::gcdist(x=c(MontBLatLon[2],locs$Lon),y=c(MontBLatLon[1],locs$Lat))
dists<-dists[1,2:(dim(dists)[2])]
MontBInd<-which(dists==min(dists))

dists<-ncf::gcdist(x=c(CarBLatLon[2],locs$Lon),y=c(CarBLatLon[1],locs$Lat))
dists<-dists[1,2:(dim(dists)[2])]
CarBInd<-which(dists==min(dists))

dists<-ncf::gcdist(x=c(MorBLatLon[2],locs$Lon),y=c(MorBLatLon[1],locs$Lat))
dists<-dists[1,2:(dim(dists)[2])]
MorBInd<-which(dists==min(dists))

dists<-ncf::gcdist(x=c(PtConcLatLon[2],locs$Lon),y=c(PtConcLatLon[1],locs$Lat))
dists<-dists[1,2:(dim(dists)[2])]
PtConcInd<-which(dists==min(dists))

dists<-ncf::gcdist(x=c(SBLatLon[2],locs$Lon),y=c(SBLatLon[1],locs$Lat))
dists<-dists[1,2:(dim(dists)[2])]
SBInd<-which(dists==min(dists))

dists<-ncf::gcdist(x=c(OxLatLon[2],locs$Lon),y=c(OxLatLon[1],locs$Lat))
dists<-dists[1,2:(dim(dists)[2])]
OxInd<-which(dists==min(dists))

dists<-ncf::gcdist(x=c(LALatLon[2],locs$Lon),y=c(LALatLon[1],locs$Lat))
dists<-dists[1,2:(dim(dists)[2])]
LAInd<-which(dists==min(dists))

#construct a region from Pt Conception to Oxnard
SBlocstouse<-PtConcInd:OxInd
plot(locs$Lon,locs$Lat,type="p",pch=20,cex=.5)
points(locs$Lon[SBlocstouse],locs$Lat[SBlocstouse],type="p",pch=20,cex=.5,col="red")  
length(SBlocstouse) #60 locations
h<-c(SBlocstouse[1],SBlocstouse[length(SBlocstouse)])
ncf::gcdist(locs$Lon[h],locs$Lat[h]) #103.8km

#First look at the whole region from Carmel Bay to Morro Bay
MorBInd-CarBInd+1 #164 locations
h<-c(CarBInd,MorBInd)
ncf::gcdist(locs$Lon[h],locs$Lat[h]) #152.8km
#it's a bit too long (compared to the above region) and has too many sites, so split it

CC1locstouse<-CarBInd:floor(mean(c(CarBInd,MorBInd)))
CC2locstouse<-ceiling(mean(c(CarBInd,MorBInd))):MorBInd
length(CC1locstouse) #82 sites
length(CC2locstouse) #82 sites
h<-c(CC1locstouse[1],CC1locstouse[length(CC1locstouse)])
ncf::gcdist(locs$Lon[h],locs$Lat[h]) #73.35km
h<-c(CC2locstouse[1],CC2locstouse[length(CC1locstouse)])
ncf::gcdist(locs$Lon[h],locs$Lat[h]) #78.46km

#So my three roughly comparable comparable regions are SBlocstouse, CC1locstouse, CC2locstouse

#***
#Pick lags for each region
#***

CC1_lag_res<-picklags(y=kelp[CC1locstouse,],x1=NO3[CC1locstouse,],x2=waves[CC1locstouse,],maxlags=c(28,8,1))
CC1_lag_res_sortAIC<-CC1_lag_res[order(CC1_lag_res$AIC),]
CC1_lag_res_sortBIC<-CC1_lag_res[order(CC1_lag_res$BIC),]
head(CC1_lag_res_sortAIC)
head(CC1_lag_res_sortBIC)
indCC1<-which(CC1_lag_res$lagAR==4 & CC1_lag_res$lag1==1 & CC1_lag_res$lag2==0)
CC1_lag_res[indCC1,]
CC1lags<-CC1_lag_res_sortBIC[1,1:3]

CC2_lag_res<-picklags(y=kelp[CC2locstouse,],x1=NO3[CC2locstouse,],x2=waves[CC2locstouse,],maxlags=c(28,8,1))
CC2_lag_res_sortAIC<-CC2_lag_res[order(CC2_lag_res$AIC),]
CC2_lag_res_sortBIC<-CC2_lag_res[order(CC2_lag_res$BIC),]
head(CC2_lag_res_sortAIC)
head(CC2_lag_res_sortBIC)
indCC2<-which(CC2_lag_res$lagAR==4 & CC2_lag_res$lag1==1 & CC2_lag_res$lag2==0)
CC2_lag_res[indCC2,]
CC2lags<-CC2_lag_res_sortBIC[1,1:3]

SB_lag_res<-picklags(y=kelp[SBlocstouse,],x1=NO3[SBlocstouse,],x2=waves[SBlocstouse,],maxlags=c(28,8,1))
SB_lag_res_sortAIC<-SB_lag_res[order(SB_lag_res$AIC),]
SB_lag_res_sortBIC<-SB_lag_res[order(SB_lag_res$BIC),]
head(SB_lag_res_sortAIC)
head(SB_lag_res_sortBIC)
indSB<-which(SB_lag_res$lagAR==4 & SB_lag_res$lag1==1 & SB_lag_res$lag2==0)
SB_lag_res[indSB,]
SBlags<-SB_lag_res_sortBIC[1,1:3]

#***
#save results
#***

save(CC1_lag_res,CC1_lag_res_sortBIC,indCC1,CC1lags,
     CC2_lag_res,CC2_lag_res_sortBIC,indCC2,CC2lags,
     SB_lag_res,SB_lag_res_sortBIC,indSB,SBlags,file=paste0(resloc,"MostResults.RData"))
