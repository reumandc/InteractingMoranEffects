#Does some linear models to confirm earlier research indicating influences of NO3 and waves 
#on kelp. Also finishes the data cleaning.
#
#The linear models just involve regressing kelp against waves and NO3 from the current time
#and last quarter, getting R^2, and then comparing to the same using surrogate data.

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): wsyn, parallel

source("MnCoefDet_Tiny2.R")
source("MnCor_Tiny.R")

#***
#Locations for storing results, and other prep
#***

resloc1<-"../Results/Kelp_DataAfterAllCleaning/"
if (!dir.exists(resloc1))
{
  dir.create(resloc1,recursive=TRUE)
}

resloc2<-"../Results/Kelp_LinearModelResults/"
if (!dir.exists(resloc2))
{
  dir.create(resloc2,recursive=TRUE)
}

theseed<-101
numsurrog<-1000
saveRDS(theseed,paste0(resloc2,"theseed"))
saveRDS(numsurrog,paste0(resloc2,"numsurrog"))

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
#Get surrogates and create centered, appropriately normalized versions of everything
#***

#centered (but not detrended) versions of the environmental variables, necessary to get surrog to work below
NO3_c<-wsyn::cleandat(NO3,times,1)$cdat
waves_c<-wsyn::cleandat(waves,times,1)$cdat

#Get surrogates for NO3 and waves - use joint surrogates of both variables, synchrony-preserving, aaft.
#We do NOT detrend the data before surrogating because it makes the surrogates less good.
set.seed(theseed) #for repeatability
h<-wsyn::surrog(dat=rbind(NO3_c,waves_c),nsurrogs=numsurrog,surrtype="aaft",syncpres=TRUE)
NO3_c_s<-lapply(FUN=function(m){return(m[1:numts,])},X=h)
waves_c_s<-lapply(FUN=function(m){return(m[(numts+1):(2*numts),])},X=h)

#Now detrend but don't variance standardize the NO3 and waves data and surrogates, and make a list for NO3 
#and one for waves with the data and then the surrogates, these are the data (and surrogates) I will then 
#work with. 
NO3<-c(list(NO3_c),NO3_c_s)
NO3<-parallel::mclapply(FUN=function(x){wsyn::cleandat(x,times=times,clev=2)$cdat},X=NO3,mc.cores=10)
waves<-c(list(waves_c),waves_c_s)
waves<-parallel::mclapply(FUN=function(x){wsyn::cleandat(x,times=times,clev=2)$cdat},X=waves,mc.cores=10)

#Option 1: Detrend and variance standardize the kelp
#kelp<-wsyn::cleandat(kelp,times,3)$cdat

#Option 2: Divide each kelp time series by its max and then detrend.
#for (counter in 1:(dim(kelp)[1]))
#{
#  kelp[counter,]<-kelp[counter,]/max(kelp[counter,])
#}
#kelp<-wsyn::cleandat(kelp,times,2)$cdat

#Option 3: Compute a kelp density per unit useable habitat. To do this, take the max kelp area value ever
#achieved for each habitat patch, use that as the measure of available habitat, and divide the kelp biomass
#time series at that location through by that value. Then detrend.
for (counter in 1:(dim(kelp)[1]))
{
  kelp[counter,]<-kelp[counter,]/(unname(quantile(kelpA[counter,],prob=.9))) #use the 90th percentile for robustness, instead of the max
}
kelp<-wsyn::cleandat(kelp,times,2)$cdat

#detrend climate indices
climinds$NPGO<-wsyn::cleandat(climinds$NPGO,times,2)$cdat
climinds$MEI<-wsyn::cleandat(climinds$MEI,times,2)$cdat
climinds$PDO<-wsyn::cleandat(climinds$PDO,times,2)$cdat

#throw out large variables you no longer need, to save memory
rm(h,NO3_c,NO3_c_s,waves_c,waves_c_s)

#***
#Save the final version of cleaned data
#***

saveRDS(kelp,file=paste0(resloc1,"Kelp_Quarterly_CleanedFinal.Rds"))
write.table(kelp,file=paste0(resloc1,"Kelp_Quarterly_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(kelpA,file=paste0(resloc1,"KelpA_Quarterly_CleanedFinal.Rds"))
write.table(kelpA,file=paste0(resloc1,"KelpA_Quarterly_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(NO3[[1]],file=paste0(resloc1,"NO3_Quarterly_CleanedFinal.Rds"))
write.table(NO3[[1]],file=paste0(resloc1,"NO3_Quarterly_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(waves[[1]],file=paste0(resloc1,"Waves_Quarterly_CleanedFinal.Rds"))
write.table(waves[[1]],file=paste0(resloc1,"Waves_Quarterly_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(locs,file=paste0(resloc1,"Locs_CleanedFinal.Rds"))
write.table(locs,file=paste0(resloc1,"Locs_CleanedFinal.csv"),row.names = FALSE,col.names = TRUE,sep=",")

saveRDS(quarters,file=paste0(resloc1,"Quarters_CleanedFinal.Rds"))
write.table(quarters,file=paste0(resloc1,"Quarters_CleanedFinal.csv"),row.names = FALSE,col.names = TRUE,sep=",")

saveRDS(climinds,file=paste0(resloc1,"Climinds_Quarterly_CleanedFinal.Rds"))
write.table(climinds,file=paste0(resloc1,"Climinds_Quarterly_CleanedFinal.csv"),row.names = FALSE,col.names = TRUE,sep=",")

#***
#Now do the linear regressions for all locations
#***

allres<-matrix(NA,length(NO3),5)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allres[counter,]<-mncoefdet_tiny2(kelp,NO3[[counter]],waves[[counter]])
}
colnames(allres)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allres,file=paste0(resloc2,"allres"))

#Thus allfsgtd gives, for each of the 5 statistics (Rsq, lag 0 NO3 coefficient, lag 1 NO3 coefficient, lag 0 waves coefficient,
#lag 1 waves coefficient), the fraction of values of that statistic, on surrogates, which were greater then the value of the 
#same statistic on data
allfsgtd<-NA*numeric(5)
for (counter in 1:5)
{
  allfsgtd[counter]<-sum(allres[2:(dim(allres)[1]),counter]>allres[1,counter])/(dim(allres)[1]-1)
}
names(allfsgtd)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allfsgtd,file=paste0(resloc2,"allfsgtd"))

#***
#Now do the same linear regressions, but for only selected geographic regions
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

#***Now do central cal, region 1

allres_CC1<-matrix(NA,length(NO3),5)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allres_CC1[counter,]<-mncoefdet_tiny2(kelp[CC1locstouse,],NO3[[counter]][CC1locstouse,],waves[[counter]][CC1locstouse,])
}
colnames(allres_CC1)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allres_CC1,file=paste0(resloc2,"allres_CC1"))

#Thus allfsgtd gives, for each of the 5 statistics (Rsq, lag 0 NO3 coefficient, lag 1 NO3 coefficient, lag 0 waves coefficient,
#lag 1 waves coefficient), the fraction of values of that statistic, on surrogates, which were greater then the value of the 
#same statistic on data
allfsgtd_CC1<-NA*numeric(5)
for (counter in 1:5)
{
  allfsgtd_CC1[counter]<-sum(allres_CC1[2:(dim(allres_CC1)[1]),counter]>allres_CC1[1,counter])/(dim(allres_CC1)[1]-1)
}
names(allfsgtd_CC1)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allfsgtd_CC1,file=paste0(resloc2,"allfsgtd_CC1"))
#allfsgtd_CC1<-readRDS(paste0(resloc2,"allfsgtd_CC1"))

#***Now do central cal, region 2

allres_CC2<-matrix(NA,length(NO3),5)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allres_CC2[counter,]<-mncoefdet_tiny2(kelp[CC2locstouse,],NO3[[counter]][CC2locstouse,],waves[[counter]][CC2locstouse,])
}
colnames(allres_CC2)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allres_CC2,file=paste0(resloc2,"allres_CC2"))

#Thus allfsgtd gives, for each of the 5 statistics (Rsq, lag 0 NO3 coefficient, lag 1 NO3 coefficient, lag 0 waves coefficient,
#lag 1 waves coefficient), the fraction of values of that statistic, on surrogates, which were greater then the value of the 
#same statistic on data
allfsgtd_CC2<-NA*numeric(5)
for (counter in 1:5)
{
  allfsgtd_CC2[counter]<-sum(allres_CC2[2:(dim(allres_CC2)[1]),counter]>allres_CC2[1,counter])/(dim(allres_CC2)[1]-1)
}
names(allfsgtd_CC2)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allfsgtd_CC2,file=paste0(resloc2,"allfsgtd_CC2"))

#***Now do the two central cal regions combined

allres_CC<-matrix(NA,length(NO3),5)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allres_CC[counter,]<-mncoefdet_tiny2(kelp[c(CC1locstouse,CC2locstouse),],
                                       NO3[[counter]][c(CC1locstouse,CC2locstouse),],
                                       waves[[counter]][c(CC1locstouse,CC2locstouse),])
}
colnames(allres_CC)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allres_CC,file=paste0(resloc2,"allres_CC"))

#Thus allfsgtd gives, for each of the 5 statistics (Rsq, lag 0 NO3 coefficient, lag 1 NO3 coefficient, lag 0 waves coefficient,
#lag 1 waves coefficient), the fraction of values of that statistic, on surrogates, which were greater then the value of the 
#same statistic on data
allfsgtd_CC<-NA*numeric(5)
for (counter in 1:5)
{
  allfsgtd_CC[counter]<-sum(allres_CC[2:(dim(allres_CC)[1]),counter]>allres_CC[1,counter])/(dim(allres_CC)[1]-1)
}
names(allfsgtd_CC)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allfsgtd_CC,file=paste0(resloc2,"allfsgtd_CC"))

#***Now do southern cal

allres_SB<-matrix(NA,length(NO3),5)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allres_SB[counter,]<-mncoefdet_tiny2(kelp[SBlocstouse,],NO3[[counter]][SBlocstouse,],waves[[counter]][SBlocstouse,])
}
colnames(allres_SB)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allres_SB,file=paste0(resloc2,"allres_SB"))

#Thus allfsgtd gives, for each of the 5 statistics (Rsq, lag 0 NO3 coefficient, lag 1 NO3 coefficient, lag 0 waves coefficient,
#lag 1 waves coefficient), the fraction of values of that statistic, on surrogates, which were greater then the value of the 
#same statistic on data
allfsgtd_SB<-NA*numeric(5)
for (counter in 1:5)
{
  allfsgtd_SB[counter]<-sum(allres_SB[2:(dim(allres_SB)[1]),counter]>allres_SB[1,counter])/(dim(allres_SB)[1]-1)
}
names(allfsgtd_SB)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allfsgtd_SB,file=paste0(resloc2,"allfsgtd_SB"))

#***Now do our three regions combined

allinds<-c(CC1locstouse,CC2locstouse,SBlocstouse)

allres_RE<-matrix(NA,length(NO3),5)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allres_RE[counter,]<-mncoefdet_tiny2(kelp[allinds,],NO3[[counter]][allinds,],waves[[counter]][allinds,])
}
colnames(allres_RE)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allres_RE,file=paste0(resloc2,"allres_RE"))

#Thus allfsgtd gives, for each of the 5 statistics (Rsq, lag 0 NO3 coefficient, lag 1 NO3 coefficient, lag 0 waves coefficient,
#lag 1 waves coefficient), the fraction of values of that statistic, on surrogates, which were greater then the value of the 
#same statistic on data
allfsgtd_RE<-NA*numeric(5)
for (counter in 1:5)
{
  allfsgtd_RE[counter]<-sum(allres_RE[2:(dim(allres_RE)[1]),counter]>allres_RE[1,counter])/(dim(allres_RE)[1]-1)
}
names(allfsgtd_RE)<-c("Rsq","NO3_l0","NO3_l1","waves_l0","waves_l1")
saveRDS(allfsgtd_RE,file=paste0(resloc2,"allfsgtd_RE"))

#***
#Now just do some lagged correlation analyses, maybe it's simpler
#***

allcorres<-matrix(NA,length(NO3),8)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allcorres[counter,]<-c(mncor_tiny(kelp,NO3[[counter]],0),
                         mncor_tiny(kelp,NO3[[counter]],1),
                         mncor_tiny(kelp,NO3[[counter]],2),
                         mncor_tiny(kelp,NO3[[counter]],3),
                         mncor_tiny(kelp,waves[[counter]],0),
                         mncor_tiny(kelp,waves[[counter]],1),
                         mncor_tiny(kelp,waves[[counter]],2),
                         mncor_tiny(kelp,waves[[counter]],3))
}
saveRDS(allcorres,paste0(resloc2,"allcorres"))

allcorfsgtd<-NA*numeric(dim(allcorres)[2])
for (counter in 1:(dim(allcorres)[2]))
{
  allcorfsgtd[counter]<-sum(allcorres[2:dim(allcorres)[1],counter]>allcorres[1,counter])/(dim(allcorres)[1]-1)
}
saveRDS(allcorfsgtd,paste0(resloc2,"allcorfsgtd"))

#now do it for just CC1

allcorres_CC1<-matrix(NA,length(NO3),8)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allcorres_CC1[counter,]<-c(mncor_tiny(kelp[CC1locstouse,],NO3[[counter]][CC1locstouse,],0),
                         mncor_tiny(kelp[CC1locstouse,],NO3[[counter]][CC1locstouse,],1),
                         mncor_tiny(kelp[CC1locstouse,],NO3[[counter]][CC1locstouse,],2),
                         mncor_tiny(kelp[CC1locstouse,],NO3[[counter]][CC1locstouse,],3),
                         mncor_tiny(kelp[CC1locstouse,],waves[[counter]][CC1locstouse,],0),
                         mncor_tiny(kelp[CC1locstouse,],waves[[counter]][CC1locstouse,],1),
                         mncor_tiny(kelp[CC1locstouse,],waves[[counter]][CC1locstouse,],2),
                         mncor_tiny(kelp[CC1locstouse,],waves[[counter]][CC1locstouse,],3))
}
saveRDS(allcorres_CC1,paste0(resloc2,"allcorres_CC1"))

allcorfsgtd_CC1<-NA*numeric(dim(allcorres_CC1)[2])
for (counter in 1:(dim(allcorres_CC1)[2]))
{
  allcorfsgtd_CC1[counter]<-sum(allcorres_CC1[2:dim(allcorres_CC1)[1],counter]>allcorres_CC1[1,counter])/(dim(allcorres_CC1)[1]-1)
}
saveRDS(allcorfsgtd_CC1,paste0(resloc2,"allcorfsgtd_CC1"))

#now do it for just CC2

allcorres_CC2<-matrix(NA,length(NO3),8)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allcorres_CC2[counter,]<-c(mncor_tiny(kelp[CC2locstouse,],NO3[[counter]][CC2locstouse,],0),
                             mncor_tiny(kelp[CC2locstouse,],NO3[[counter]][CC2locstouse,],1),
                             mncor_tiny(kelp[CC2locstouse,],NO3[[counter]][CC2locstouse,],2),
                             mncor_tiny(kelp[CC2locstouse,],NO3[[counter]][CC2locstouse,],3),
                             mncor_tiny(kelp[CC2locstouse,],waves[[counter]][CC2locstouse,],0),
                             mncor_tiny(kelp[CC2locstouse,],waves[[counter]][CC2locstouse,],1),
                             mncor_tiny(kelp[CC2locstouse,],waves[[counter]][CC2locstouse,],2),
                             mncor_tiny(kelp[CC2locstouse,],waves[[counter]][CC2locstouse,],3))
}
saveRDS(allcorres_CC2,paste0(resloc2,"allcorres_CC2"))

allcorfsgtd_CC2<-NA*numeric(dim(allcorres_CC2)[2])
for (counter in 1:(dim(allcorres_CC2)[2]))
{
  allcorfsgtd_CC2[counter]<-sum(allcorres_CC2[2:dim(allcorres_CC2)[1],counter]>allcorres_CC2[1,counter])/(dim(allcorres_CC2)[1]-1)
}
saveRDS(allcorfsgtd_CC2,paste0(resloc2,"allcorfsgtd_CC2"))

#now do it for CC only (CC1+CC2)

inds<-c(CC1locstouse,CC2locstouse)
allcorres_CC<-matrix(NA,length(NO3),8)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allcorres_CC[counter,]<-c(mncor_tiny(kelp[inds,],NO3[[counter]][inds,],0),
                             mncor_tiny(kelp[inds,],NO3[[counter]][inds,],1),
                             mncor_tiny(kelp[inds,],NO3[[counter]][inds,],2),
                             mncor_tiny(kelp[inds,],NO3[[counter]][inds,],3),
                             mncor_tiny(kelp[inds,],waves[[counter]][inds,],0),
                             mncor_tiny(kelp[inds,],waves[[counter]][inds,],1),
                             mncor_tiny(kelp[inds,],waves[[counter]][inds,],2),
                             mncor_tiny(kelp[inds,],waves[[counter]][inds,],3))
}
saveRDS(allcorres_CC,paste0(resloc2,"allcorres_CC"))

allcorfsgtd_CC<-NA*numeric(dim(allcorres_CC)[2])
for (counter in 1:(dim(allcorres_CC)[2]))
{
  allcorfsgtd_CC[counter]<-sum(allcorres_CC[2:dim(allcorres_CC)[1],counter]>allcorres_CC[1,counter])/(dim(allcorres_CC)[1]-1)
}
saveRDS(allcorfsgtd_CC,paste0(resloc2,"allcorfsgtd_CC"))

#now do it for just the SB region

allcorres_SB<-matrix(NA,length(NO3),8)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allcorres_SB[counter,]<-c(mncor_tiny(kelp[SBlocstouse,],NO3[[counter]][SBlocstouse,],0),
                             mncor_tiny(kelp[SBlocstouse,],NO3[[counter]][SBlocstouse,],1),
                            mncor_tiny(kelp[SBlocstouse,],NO3[[counter]][SBlocstouse,],2),
                            mncor_tiny(kelp[SBlocstouse,],NO3[[counter]][SBlocstouse,],3),
                             mncor_tiny(kelp[SBlocstouse,],waves[[counter]][SBlocstouse,],0),
                             mncor_tiny(kelp[SBlocstouse,],waves[[counter]][SBlocstouse,],1),
                            mncor_tiny(kelp[SBlocstouse,],waves[[counter]][SBlocstouse,],2),
                            mncor_tiny(kelp[SBlocstouse,],waves[[counter]][SBlocstouse,],3))}
saveRDS(allcorres_SB,paste0(resloc2,"allcorres_SB"))

allcorfsgtd_SB<-NA*numeric(dim(allcorres_SB)[2])
for (counter in 1:(dim(allcorres_SB)[2]))
{
  allcorfsgtd_SB[counter]<-sum(allcorres_SB[2:dim(allcorres_SB)[1],counter]>allcorres_SB[1,counter])/(dim(allcorres_SB)[1]-1)
}
saveRDS(allcorfsgtd_SB,paste0(resloc2,"allcorfsgtd_SB"))

#now do it for the combination of the three regions

inds<-c(CC1locstouse,CC2locstouse,SBlocstouse)
allcorres_RE<-matrix(NA,length(NO3),8)
for (counter in 1:length(NO3))
{
  print(paste(counter,"of",length(NO3)))
  allcorres_RE[counter,]<-c(mncor_tiny(kelp[inds,],NO3[[counter]][inds,],0),
                            mncor_tiny(kelp[inds,],NO3[[counter]][inds,],1),
                            mncor_tiny(kelp[inds,],NO3[[counter]][inds,],2),
                            mncor_tiny(kelp[inds,],NO3[[counter]][inds,],3),
                            mncor_tiny(kelp[inds,],waves[[counter]][inds,],0),
                            mncor_tiny(kelp[inds,],waves[[counter]][inds,],1),
                            mncor_tiny(kelp[inds,],waves[[counter]][inds,],2),
                            mncor_tiny(kelp[inds,],waves[[counter]][inds,],3))}
saveRDS(allcorres_RE,paste0(resloc2,"allcorres_RE"))

allcorfsgtd_RE<-NA*numeric(dim(allcorres_RE)[2])
for (counter in 1:(dim(allcorres_RE)[2]))
{
  allcorfsgtd_RE[counter]<-sum(allcorres_RE[2:dim(allcorres_RE)[1],counter]>allcorres_RE[1,counter])/(dim(allcorres_RE)[1]-1)
}
saveRDS(allcorfsgtd_RE,paste0(resloc2,"allcorfsgtd_RE"))

#On balance, what all this justifies is probably using lags 0 and 1 for NO3 and lag 0 only for waves. 
#Probably best to present just the lagged correlation analyses, for simplicity, and because we don't really learn anything 
#from the regression results that we don't get from the correlation analyses, namely, what lags to use for NO3 and waves.
#
#Just justify a lag of 4 or 8 or 12 for kelp on the grounds that it's a perrenial. Do subsequent analyses with all three
#of these kelp lags, to test for sensitivity of results to that.