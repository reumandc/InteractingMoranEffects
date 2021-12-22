#Cleans the data in a very basic way. Just removes quarters and sites with no data for one 
#of the variables, and puts sites in "paddle order," the order you would pass the sites in 
#if you paddled a kayak from north to south, going into all the bays and such.

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): profvis

#***
#Location for storing results and other prep
#***

resloc<-"../Results/Kelp_DataAfterSuperBasicCleaning/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***
#Load the data
#***

datloc<-"../Data/"

#Load kelp data and change NaNs to NAs. Use the biomass data, sites ordered as in the 
#original ordering (which should match with the environmental variables), not the paddling
#order (not yet, anyway - we reorder sites later)
kelp_df<-read.csv(file=paste0(datloc,"Kelp_Bio_2019_v3.csv"),header=TRUE)
for (counter in 1:(dim(kelp_df)[2]))
{
  h<-kelp_df[,counter]
  h[is.nan(h)]<-NA
  kelp_df[,counter]<-h
}

#make it a matrix, locations by quarters
kelp_mat<-as.matrix(kelp_df[,2:(dim(kelp_df)[2])])

#check for the same values
h1<-kelp_df[,2:(dim(kelp_df)[2])]
h2<-kelp_mat
dim(h1)
dim(h2) #dimensions should be the same
inds_df<-which(is.na(h1))
inds_mat<-which(is.na(h2))
length(inds_df)
length(inds_mat)
sum(inds_df==inds_mat) #both the matrix and data frame versions should have the same NAs
inds_neq<-which(h1!=h2)
length(inds_neq) #should be none
prod(dim(h1))
sum(h1==h2,na.rm=TRUE)+length(inds_df) #number of equal locations plus number of NA locations should be all the locations

kelp<-unname(kelp_mat) #this is a matrix consisting of numbers and NAs, dimensions are locations by quarters
rm(kelp_df,kelp_mat,h,h1,h2,inds_df,inds_mat,inds_neq,counter)

#Same for kelp area - this is area of kelp canopy
kelp_df<-read.csv(file=paste0(datloc,"Kelp_Area_2019_v3.csv"),header=TRUE)
for (counter in 1:(dim(kelp_df)[2]))
{
  h<-kelp_df[,counter]
  h[is.nan(h)]<-NA
  kelp_df[,counter]<-h
}

#make it a matrix, locations by quarters
kelp_mat<-as.matrix(kelp_df[,2:(dim(kelp_df)[2])])

#check for the same values
h1<-kelp_df[,2:(dim(kelp_df)[2])]
h2<-kelp_mat
dim(h1)
dim(h2) #dimensions should be the same
inds_df<-which(is.na(h1))
inds_mat<-which(is.na(h2))
length(inds_df)
length(inds_mat)
sum(inds_df==inds_mat) #both the matrix and data frame versions should have the same NAs
inds_neq<-which(h1!=h2)
length(inds_neq) #should be none
prod(dim(h1))
sum(h1==h2,na.rm=TRUE)+length(inds_df) #number of equal locations plus number of NA locations should be all the locations

kelpA<-unname(kelp_mat) #this is a matrix consisting of numbers and NAs, dimensions are locations by quarters
rm(kelp_df,kelp_mat,h,h1,h2,inds_df,inds_mat,inds_neq,counter)

#load the quarters
quarters<-read.csv(file=paste0(datloc,"Quarters.csv"),header=TRUE)
head(quarters)

#check for dimensional consistency with kelp
dim(quarters)
dim(kelp)

#load the locations
sites<-read.csv(file=paste0(datloc,"Sites.csv"),header=TRUE)

#check for dimensional consistency with kelp
dim(sites)
dim(kelp)

#load nitrates
NO3<-read.csv(file=paste0(datloc,"NO3_Mean_2019.csv"),header=TRUE)
NO3<-as.matrix(NO3[,2:(dim(NO3)[2])])
NO3[is.nan(NO3)]<-NA
dim(NO3)
sum(!is.finite(NO3))
sum(is.na(NO3))
NO3<-unname(NO3)

#load max wave height
wave<-read.csv(file=paste0(datloc,"Hs_Max_2019.csv"),head=TRUE)
wave<-as.matrix(wave[,2:(dim(wave)[2])])
wave[is.nan(wave)]<-NA
dim(wave)
sum(!is.finite(wave))
sum(is.na(wave))
wave<-unname(wave)

#paddle ordering - encodes an ordering of the sites which is the order you would pass them in 
#if you paddled north to south along the coast in a kayak, going into all the bays and such
paddle<-read.csv(file=paste0(datloc,"Paddle_coords.csv"),head=TRUE)
paddle$Paddle_Number<-1:dim(paddle)[1]
paddle<-paddle[order(paddle$Site_Number),]

#check for consistency with sites
sum(paddle$Lat==sites$Lat)
sum(paddle$Lon==sites$Lon)
rm(sites)

#climatic indices
climind<-read.csv(file=paste0(datloc,"Indices.csv"),head=TRUE)
head(climind)
dim(climind)
dim(kelp)

#***
#put everything in paddle order and remove redundancies
#***

paddle<-paddle[order(paddle$Paddle_Number),]

kelp<-kelp[paddle$Site_Number,]
kelpA<-kelpA[paddle$Site_Number,]
NO3<-NO3[paddle$Site_Number,]
wave<-wave[paddle$Site_Number,]
locs<-paddle[,1:2]
rm(paddle)

#***
#check visually to make sure it is paddle order
#***

pl<-.1
plot(locs$Lon[1],locs$Lat[1],type="p",pch=20,cex=.1,ylim=range(locs$Lat),xlim=range(locs$Lon))
profvis::pause(pl)
for (counter in 2:(dim(locs)[1]))
{
  points(locs$Lon[counter],locs$Lat[counter],type="p",pch=20,cex=.1)
  profvis::pause(pl)
}
#paddling order visually confirmed

#***
#Do some early-stage cleaning
#***

#the first 12 quarters were all NAs for the waves, so throw them out for everything
kelp<-kelp[,-(1:12)]
kelpA<-kelpA[,-(1:12)]
NO3<-NO3[,-(1:12)]
wave<-wave[,-(1:12)]
quarters<-quarters[-(1:12),]
climind<-climind[-(1:12),]

#look for other quarters we should maybe omit
apply(FUN=function(x){sum(is.na(x))},X=kelp,MARGIN=2)
apply(FUN=function(x){sum(is.na(x))},X=NO3,MARGIN=2)
apply(FUN=function(x){sum(is.na(x))},X=wave,MARGIN=2)
#decided to keep all available data for all remaining quarters

#same for sites
hN<-apply(FUN=function(x){sum(is.na(x))},X=NO3,MARGIN=1)
hN
hW<-apply(FUN=function(x){sum(is.na(x))},X=wave,MARGIN=1)
hW
#these aren't too bad

#what about kelp, how many NAs, how many zeros?
numNAk<-apply(FUN=function(x){sum(is.na(x))},X=kelp,MARGIN=1)
numNAk
num0k<-apply(FUN=function(x){sum(x==0,na.rm=TRUE)},X=kelp,MARGIN=1)
num0k

#some sites have no kelp data, so remove those sites for all variables
kelp<-kelp[numNAk!=132,]
kelpA<-kelpA[numNAk!=132,]
NO3<-NO3[numNAk!=132,]
wave<-wave[numNAk!=132,]
locs<-locs[numNAk!=132,]
rm(hN,hW)

#***
#one more quick check of basic dimensions
#***

dim(kelp)
dim(kelpA)
dim(NO3)
dim(wave)
dim(locs)
dim(quarters)
dim(climind)

#***
#save the results
#***

saveRDS(kelp,file=paste0(resloc,"Kelp_Quarterly_CleanedSuperBasic.Rds"))
write.table(kelp,file=paste0(resloc,"Kelp_Quarterly_CleanedSuperBasic.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(kelpA,file=paste0(resloc,"KelpA_Quarterly_CleanedSuperBasic.Rds"))
write.table(kelpA,file=paste0(resloc,"KelpA_Quarterly_CleanedSuperBasic.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(NO3,file=paste0(resloc,"NO3_Quarterly_CleanedSuperBasic.Rds"))
write.table(NO3,file=paste0(resloc,"NO3_Quarterly_CleanedSuperBasic.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(wave,file=paste0(resloc,"Waves_Quarterly_CleanedSuperBasic.Rds"))
write.table(wave,file=paste0(resloc,"Waves_Quarterly_CleanedSuperBasic.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(locs,file=paste0(resloc,"Locs_CleanedSuperBasic.Rds"))
write.table(locs,file=paste0(resloc,"Locs_CleanedSuperBasic.csv"),row.names = FALSE,col.names = TRUE,sep=",")

saveRDS(quarters,file=paste0(resloc,"Quarters_CleanedSuperBasic.Rds"))
write.table(quarters,file=paste0(resloc,"Quarters_CleanedSuperBasic.csv"),row.names = FALSE,col.names = TRUE,sep=",")

saveRDS(climind,file=paste0(resloc,"Climinds_Quarterly_CleanedSuperBasic.Rds"))
write.table(climind,file=paste0(resloc,"Climinds_Quarterly_CleanedSuperBasic.csv"),row.names = FALSE,col.names = TRUE,sep=",")


