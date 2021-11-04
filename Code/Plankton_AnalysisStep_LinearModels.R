#Does some correlation analyses to confirm influences of temp and calfin on pci, and to get lags. 
#Also finishes the data cleaning. The linear models done for kelp are skipped here, because we
#ended up not needing them for kelp.

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): wsyn, parallel

source("MnCor_Tiny.R")

#***
#Locations for storing results, and other prep
#***

resloc1<-"../Results/Plankton_DataAfterAllCleaning/"
if (!dir.exists(resloc1))
{
  dir.create(resloc1,recursive=TRUE)
}

resloc2<-"../Results/Plankton_LinearModelResults/"
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

datloc<-"../Results/Plankton_DataAfterBasicCleaning/"
pci<-readRDS(paste0(datloc,"Pci_CleanedBasic.Rds")) 
calfin<-readRDS(paste0(datloc,"Calfin_CleanedBasic.Rds")) 
temp<-readRDS(paste0(datloc,"Temp_CleanedBasic.Rds")) 
year<-readRDS(paste0(datloc,"Year_CleanedBasic.Rds"))
site<-readRDS(paste0(datloc,"Site_CleanedBasic.Rds"))

numts<-dim(pci)[1]
lents<-dim(pci)[2]

# #***
# #Get surrogates and create centered, appropriately normalized versions of everything - approach 1
# #***
# 
# #centered (but not detrended) versions of the "environmental" variables (considered temp and calfin), 
# #necessary to get surrog to work below
# calfin_c<-wsyn::cleandat(calfin,year,1)$cdat
# temp_c<-wsyn::cleandat(temp,year,1)$cdat
# 
# #Get surrogates for calfin and temp - use joint surrogates of both variables, synchrony-preserving, aaft.
# #We do NOT detrend the data before surrogating because it makes the surrogates less good.
# set.seed(theseed) #for repeatability
# h<-wsyn::surrog(dat=rbind(calfin_c,temp_c),nsurrogs=numsurrog,surrtype="aaft",syncpres=TRUE)
# calfin_c_s<-lapply(FUN=function(m){return(m[1:numts,])},X=h)
# temp_c_s<-lapply(FUN=function(m){return(m[(numts+1):(2*numts),])},X=h)
# 
# #Now detrend but don't variance standardize the calfin and temp data and surrogates, and make a list for calfin 
# #and one for temp with the data and then the surrogates, these are the data (and surrogates) I will then 
# #work with. 
# calfin<-c(list(calfin_c),calfin_c_s)
# calfin<-parallel::mclapply(FUN=function(x){wsyn::cleandat(x,times=year,clev=2)$cdat},X=calfin,mc.cores=10)
# temp<-c(list(temp_c),temp_c_s)
# temp<-parallel::mclapply(FUN=function(x){wsyn::cleandat(x,times=year,clev=2)$cdat},X=temp,mc.cores=10)
# 
# #detrend but don't variance standardize the pci data
# pci<-wsyn::cleandat(pci,year,2)$cdat
# 
# #throw out large variables you no longer need, to save memory
# rm(h,calfin_c,calfin_c_s,temp_c,temp_c_s)

#***
#Cleaning approach 2
#***

calfin_cd<-wsyn::cleandat(calfin,times=1:lents,clev=2)
calfin<-calfin_cd$cdat
pci<-wsyn::cleandat(pci,times=1:lents,clev=2)$cdat
temp<-wsyn::cleandat(temp,times=1:lents,clev=2)$cdat

saveRDS(pci,file=paste0(resloc1,"Pci_CleanedFinal.Rds"))
write.table(pci,file=paste0(resloc1,"Pci_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(calfin,file=paste0(resloc1,"Calfin_CleanedFinal.Rds"))
write.table(calfin,file=paste0(resloc1,"Calfin_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(temp,file=paste0(resloc1,"Temp_CleanedFinal.Rds"))
write.table(temp,file=paste0(resloc1,"Temp_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(year,file=paste0(resloc1,"Year_CleanedFinal.Rds"))
write.table(year,file=paste0(resloc1,"Year_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(site,file=paste0(resloc1,"Site_CleanedFinal.Rds"))
write.table(site,file=paste0(resloc1,"Site_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

#***
#Save the final version of cleaned data
#***

saveRDS(pci,file=paste0(resloc1,"Pci_CleanedFinal.Rds"))
write.table(pci,file=paste0(resloc1,"Pci_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(calfin[[1]],file=paste0(resloc1,"Calfin_CleanedFinal.Rds"))
write.table(calfin[[1]],file=paste0(resloc1,"Calfin_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(temp[[1]],file=paste0(resloc1,"Temp_CleanedFinal.Rds"))
write.table(temp[[1]],file=paste0(resloc1,"Temp_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(year,file=paste0(resloc1,"Year_CleanedFinal.Rds"))
write.table(year,file=paste0(resloc1,"Year_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

saveRDS(site,file=paste0(resloc1,"Site_CleanedFinal.Rds"))
write.table(site,file=paste0(resloc1,"Site_CleanedFinal.csv"),row.names = FALSE,col.names = FALSE,sep=",")

#***
#Now do some lagged correlation analyses
#***

#None of the below is getting where I want to go at the moment

allcorres<-matrix(NA,length(calfin),20)
for (counter in 1:length(calfin))
{
  print(paste(counter,"of",length(calfin)))
  allcorres[counter,]<-c(mncor_tiny(pci,calfin[[counter]],0),
                         mncor_tiny(pci,calfin[[counter]],1),
                         mncor_tiny(pci,calfin[[counter]],2),
                         mncor_tiny(pci,calfin[[counter]],3),
                         mncor_tiny(pci,calfin[[counter]],4),
                         mncor_tiny(pci,calfin[[counter]],5),
                         mncor_tiny(pci,calfin[[counter]],6),
                         mncor_tiny(pci,calfin[[counter]],7),
                         mncor_tiny(pci,calfin[[counter]],8),
                         mncor_tiny(pci,calfin[[counter]],9),
                         mncor_tiny(pci,temp[[counter]],0),
                         mncor_tiny(pci,temp[[counter]],1),
                         mncor_tiny(pci,temp[[counter]],2),
                         mncor_tiny(pci,temp[[counter]],3),
                         mncor_tiny(pci,temp[[counter]],4),
                         mncor_tiny(pci,temp[[counter]],5),
                         mncor_tiny(pci,temp[[counter]],6),
                         mncor_tiny(pci,temp[[counter]],7),
                         mncor_tiny(pci,temp[[counter]],8),
                         mncor_tiny(pci,temp[[counter]],9))
}
saveRDS(allcorres,paste0(resloc2,"allcorres"))
#allcorres<-readRDS(file=paste0(resloc2,"allcorres"))

allcorfsgtd<-NA*numeric(dim(allcorres)[2])
for (counter in 1:(dim(allcorres)[2]))
{
  allcorfsgtd[counter]<-sum(allcorres[2:dim(allcorres)[1],counter]>allcorres[1,counter])/(dim(allcorres)[1]-1)
}
saveRDS(allcorfsgtd,paste0(resloc2,"allcorfsgtd"))
#allcorfsgtd<-readRDS(file=paste0(resloc2,"allcorfsgtd"))

#OK, this illuminates what kind of lag structure I might use for calfin, but not for temp.
#Recall, there was a quarter-phase lag between temp and pci across all long (>4yr) timescales,
#and an antiphase relationship between cal fin and pci. The anti-phase relationship seems
#to translate sensibly to a lag-0 dependence of pci on calfin with a negative coefficient.
#But not sure what to do to translate the wavelet coherence result between temp and pci into
#an ARMA-type model.

#To get started thinking, let's try out the same regression approach using 0 and 1 lags as 
#was used for kelp

allres<-matrix(NA,length(calfin),5)
for (counter in 1:length(calfin))
{
  print(paste(counter,"of",length(calfin)))
  allres[counter,]<-mncoefdet_tiny2(pci,calfin[[counter]],temp[[counter]])
}
colnames(allres)<-c("Rsq","calfin_l0","calfin_l1","temp_l0","temp_l1")
#saveRDS(allres,file=paste0(resloc2,"allres"))

#Thus allfsgtd gives, for each of the 5 statistics (Rsq, lag 0 calfin coefficient, lag 1 calfin 
#coefficient, lag 0 temp coefficient, lag 1 temp coefficient), the fraction of values of that 
#statistic, on surrogates, which were greater then the value of the same statistic on data
allfsgtd<-NA*numeric(5)
for (counter in 1:5)
{
  allfsgtd[counter]<-sum(allres[2:(dim(allres)[1]),counter]>allres[1,counter])/(dim(allres)[1]-1)
}
names(allfsgtd)<-c("Rsq","calfin_l0","calfin_l1","temp_l0","temp_l1")
#saveRDS(allfsgtd_CC1,file=paste0(resloc2,"allfsgtd_CC1"))

mncoefdet_tiny_experiment<-function(y,x1,x2)
{
  if (any(dim(y)!=dim(x1)) || any(dim(y)!=dim(x2)))
  {
    stop("Error in mncoefdet_tiny_experiment: y, x1 and x2 input matrices must be the same dimensions")
  }
  
  numts<-dim(y)[1]
  lents<-dim(y)[2]
  
  res_Rsq<-NA*numeric(numts)
  res_x1coefs<-matrix(NA,numts,2)
  res_x2coefs<-matrix(NA,numts,2)
  for (counter in 1:numts)
  {
    #extract time series for this location
    tsy<-y[counter,]
    tsx1<-x1[counter,]
    tsx2<-x2[counter,]
    
    #implement the desired lags
    tsy<-tsy[2:lents]
    tsx1_l0<-tsx1[2:lents]
    tsx1_l1<-tsx1[1:(lents-1)]
    tsx2_l0<-tsx2[2:lents]
    tsx2_l1<-tsx2[1:(lents-1)]
    
    #do the regressiom
    mod<-lm(tsy~tsx1_l0+tsx2_l0+tsx2_l1)
    res_Rsq[counter]<-summary(mod)$r.squared
    res_x1coefs[counter,]<-unname(coef(mod)[2:3])
    res_x2coefs[counter,]<-unname(coef(mod)[4:5])
  }
  
  return(c(mean(res_Rsq),apply(FUN=mean,X=res_NO3coefs,MARGIN=2),apply(FUN=mean,X=res_wavescoefs,MARGIN=2)))
}