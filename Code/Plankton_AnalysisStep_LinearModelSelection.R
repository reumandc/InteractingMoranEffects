#Does some linear models and information criterion analyses to select numbers of lags to use in 
#ARMA-type models of plankton.

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): ncf

source("PickLags.R")

#***
#Locations for storing results, and other prep
#***

resloc<-"../Results/Plankton_LinearModelSelection/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***
#Load the data
#***

#load the data
datloc<-"../Results/Plankton_DataAfterAllCleaning/"
pci<-readRDS(paste0(datloc,"Pci_CleanedFinal.Rds")) 
calfin<-readRDS(paste0(datloc,"Calfin_CleanedFinal.Rds")) 
temp<-readRDS(paste0(datloc,"Temp_CleanedFinal.Rds")) 
locs<-readRDS(paste0(datloc,"Site_CleanedFinal.Rds"))
year<-readRDS(paste0(datloc,"Year_CleanedFinal.Rds"))
times<-1:length(year)

numts<-dim(pci)[1]
lents<-dim(pci)[2]

#***
#Pick lags 
#***

lag_res<-picklags_mat(y=pci,x1=temp,x2=calfin,maxlags=c(1,1,2))

lag_res$AICwt<-exp(-.5*lag_res$DeltaAIC)
lag_res$AICwt<-lag_res$AICwt/sum(lag_res$AICwt)
lag_res$BICwt<-exp(-.5*lag_res$DeltaBIC)
lag_res$BICwt<-lag_res$BICwt/sum(lag_res$BICwt)

lag_res_sortAIC<-lag_res[order(lag_res$AIC),]
lag_res_sortBIC<-lag_res[order(lag_res$BIC),]
lag_res_sortAIC$AICwt_cumsum<-cumsum(lag_res_sortAIC$AICwt)
lag_res_sortBIC$BICwt_cumsum<-cumsum(lag_res_sortBIC$BICwt)

head(lag_res_sortAIC,20)
head(lag_res_sortBIC,20)

range(lag_res$lagAR[lag_res$DeltaAIC<=3])
range(lag_res$lag1[lag_res$DeltaAIC<=3])
range(lag_res$lag2[lag_res$DeltaAIC<=3])

range(lag_res$lagAR[lag_res$DeltaBIC<=3])
range(lag_res$lag1[lag_res$DeltaBIC<=3])
range(lag_res$lag2[lag_res$DeltaBIC<=3])

range(lag_res$lagAR[lag_res$DeltaAIC<=2])
range(lag_res$lag1[lag_res$DeltaAIC<=2])
range(lag_res$lag2[lag_res$DeltaAIC<=2])

range(lag_res$lagAR[lag_res$DeltaBIC<=2])
range(lag_res$lag1[lag_res$DeltaBIC<=2])
range(lag_res$lag2[lag_res$DeltaBIC<=2])
