#Does a more involved data cleaning, starting from the super basic cleaning done in 
#AnalysisStep_CleanData_SuperBasic.R. Just throws out sites with "too many" NAs/zeros,
#and fills remaining NAs with seasonal medians.

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): 

source("FillSeasonalMedian.R")
source("FindZeroNAYears.R")

#***
#Location for storing results and other prep
#***

resloc<-"../Results/DataAfterMoreInvolvedCleaning/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***
#Load the quarterly data with only super basic cleaning done to it
#***

#load the data
datloc<-"../Results/DataAfterSuperBasicCleaning/"
kelpQ<-readRDS(paste0(datloc,"Kelp_Quarterly_CleanedSuperBasic.Rds")) 
kelpAQ<-readRDS(paste0(datloc,"KelpA_Quarterly_CleanedSuperBasic.Rds")) 
NO3Q<-readRDS(paste0(datloc,"NO3_Quarterly_CleanedSuperBasic.Rds")) 
wavesQ<-readRDS(paste0(datloc,"Waves_Quarterly_CleanedSuperBasic.Rds")) 
locs<-readRDS(paste0(datloc,"Locs_CleanedSuperBasic.Rds"))
quarters<-readRDS(paste0(datloc,"Quarters_CleanedSuperBasic.Rds"))
climinds<-readRDS(paste0(datloc,"Climinds_Quarterly_CleanedSuperBasic.Rds"))
timesQ<-1:(dim(quarters)[1])

#***
#Slightly more extensive data prep - throw out sites with too many NAs or zeroes, then fill remaining NAs 
#with seasonal medians
#***

#Remove all sites which have kelp 0/NA for more than 3 whole calendar years any time during the data 
#duration. This is a conservative definition of "persistent" kelp sites - we want to work with 
#sites driven by kelp dynamics, not sediment dynamics.
badyrs<-apply(FUN=find_zero_NA_years,X=kelpQ,MARGIN=1)
hist(badyrs)
numbadyrslte<-c()
for (counter in 0:20)
{
  numbadyrslte[counter+1]<-sum(badyrs<=counter)
}
sum(badyrs==0)
sum(badyrs<=1)
sum(badyrs<=2)
sum(badyrs<=3)
plot(0:20,numbadyrslte,ylim=c(0,dim(kelpQ)[1]))
numbadyrslte<-data.frame(numlocs=0:20,numbadyrs=numbadyrslte)

goodlocs<-which(badyrs<=3)
kelpQ<-kelpQ[goodlocs,]
kelpAQ<-kelpAQ[goodlocs,]
NO3Q<-NO3Q[goodlocs,]
wavesQ<-wavesQ[goodlocs,]
locs<-locs[goodlocs,]
rownames(locs)<-1:dim(locs)[1]

#Fill remaining NAs with seasonal medians
kelp_mf<-t(apply(FUN=fill_with_seasonal_median,X=kelpQ,MARGIN=1))
sum(kelpQ==kelp_mf,na.rm=TRUE)+sum(is.na(kelpQ))
prod(dim(kelpQ))
kelpQ<-kelp_mf
rm(kelp_mf)

#do the same for kelpA
kelp_mf<-t(apply(FUN=fill_with_seasonal_median,X=kelpAQ,MARGIN=1))
sum(kelpAQ==kelp_mf,na.rm=TRUE)+sum(is.na(kelpAQ))
prod(dim(kelpAQ))
kelpAQ<-kelp_mf
rm(kelp_mf)

#do the same for nitrates and waves 
NO3Q<-t(apply(FUN=fill_with_seasonal_median,X=NO3Q,MARGIN=1))
wavesQ<-t(apply(FUN=fill_with_seasonal_median,X=wavesQ,MARGIN=1))

#check to make sure there are now no missing values
apply(FUN=function(x){sum(is.na(x))},X=kelpQ,MARGIN=2)
apply(FUN=function(x){sum(is.na(x))},X=NO3Q,MARGIN=2)
apply(FUN=function(x){sum(is.na(x))},X=wavesQ,MARGIN=2)

numlocs<-dim(kelpQ)[1]
lents<-dim(kelpQ)[2]

rm(badyrs,counter,fill_with_seasonal_median,find_zero_NA_years,goodlocs,numbadyrslte)

#***
#Save the result
#***

saveRDS(kelpQ,file=paste0(resloc,"Kelp_Quarterly_CleanedMoreInvolved.Rds"))
write.table(kelpQ,file=paste0(resloc,"Kelp_Quarterly_CleanedMoreInvolved.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(kelpAQ,file=paste0(resloc,"KelpA_Quarterly_CleanedMoreInvolved.Rds"))
write.table(kelpAQ,file=paste0(resloc,"KelpA_Quarterly_CleanedMoreInvolved.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(NO3Q,file=paste0(resloc,"NO3_Quarterly_CleanedMoreInvolved.Rds"))
write.table(NO3Q,file=paste0(resloc,"NO3_Quarterly_CleanedMoreInvolved.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(wavesQ,file=paste0(resloc,"Waves_Quarterly_CleanedMoreInvolved.Rds"))
write.table(wavesQ,file=paste0(resloc,"Waves_Quarterly_CleanedMoreInvolved.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(locs,file=paste0(resloc,"Locs_CleanedMoreInvolved.Rds"))
write.table(locs,file=paste0(resloc,"Locs_CleanedMoreInvolved.csv"),row.names = FALSE,col.names = TRUE,sep=",")

saveRDS(quarters,file=paste0(resloc,"Quarters_CleanedMoreInvolved.Rds"))
write.table(quarters,file=paste0(resloc,"Quarters_CleanedMoreInvolved.csv"),row.names = FALSE,col.names = TRUE,sep=",")

saveRDS(climinds,file=paste0(resloc,"Climinds_Quarterly_CleanedMoreInvolved.Rds"))
write.table(climinds,file=paste0(resloc,"Climinds_Quarterly_CleanedMoreInvolved.csv"),row.names = FALSE,col.names = TRUE,sep=",")





