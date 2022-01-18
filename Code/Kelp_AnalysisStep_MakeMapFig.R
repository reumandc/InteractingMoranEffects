#Makes a map figure for the paper that orients the reader to the kelp data we used

#***
#External codes needed
#***

#packages needed (invoked with "::"): 

#***
#Location for storing results and other prep
#***

resloc<-"../Results/Kelp_MapFig/"
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
#Which site is closest to pt conception, same for SB and some other locations, for later use
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

#***
#Select some regions to focus on
#***

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

#So my three roughly comparable regions are SBlocstouse, CC1locstouse, CC2locstouse. Save them.
saveRDS(SBlocstouse,file=paste0(resloc,"SBlocstouse.Rds"))
saveRDS(CC1locstouse,file=paste0(resloc,"CC1locstouse.Rds"))
saveRDS(CC2locstouse,file=paste0(resloc,"CC2locstouse.Rds"))

#***
#Organize the data into nice format for Max
#***

CC1kelpdat<-kelp[CC1locstouse,] #each row is a location with time running left to right, columns corresponding to quarters as in the variable "quarters"
CC1locs<-locs[CC1locstouse,] #lat lon of locations corresponding to rows of CC1kelpdat

CC2kelpdat<-kelp[CC2locstouse,] #each row is a location with time running left to right, columns corresponding to quarters as in the variable "quarters"
CC2locs<-locs[CC2locstouse,] #lat lon of locations corresponding to rows of CC2kelpdat

SBkelpdat<-kelp[SBlocstouse,] #each row is a location with time running left to right, columns corresponding to quarters as in the variable "quarters"
SBlocs<-locs[SBlocstouse,] #lat lon of locations corresponding to rows of SBkelpdat

#temp save files just so I can send these all to Max without him even accessing the repo
saveRDS(CC1kelpdat,file=paste0(resloc,"CC1kelpdat.Rds"))
saveRDS(CC2kelpdat,file=paste0(resloc,"CC2kelpdat.Rds"))
saveRDS(SBkelpdat,file=paste0(resloc,"SBkelpdat.Rds"))
saveRDS(CC1locs,file=paste0(resloc,"CC1locs.Rds"))
saveRDS(CC2locs,file=paste0(resloc,"CC2locs.Rds"))
saveRDS(SBlocs,file=paste0(resloc,"SBlocs.Rds"))

#OK, MAX, TAKE IT FROM HERE...