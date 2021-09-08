#Clean the plankton data. Since data were already cleaned by our earlier efforts, 
#published in Plos Comp Biol, not too much is needed here.

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): 

#***
#Location for storing results and other prep
#***

resloc<-"../Results/Plankton_DataAfterCleaning/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***
#Load the data
#***

datloc<-"../Data/PlanktonData/"

#load the phytoplankton color index (pci) data
pci<-read.csv(file=paste0(datloc,"PCIannual.csv"),header=FALSE)
pci_mat<-as.matrix(pci)
sum(pci_mat==pci)
prod(dim(pci))
pci<-pci_mat
rm(pci_mat)
year<-unname(pci[1,2:(dim(pci)[2])])
site<-pci[2:(dim(pci)[1]),1]
pci<-pci[2:(dim(pci)[1]),2:(dim(pci)[2])]

#load the cal fin dat
calfin<-read.csv(file=paste0(datloc,"Calanus_finmarchicus.csv"),header=FALSE)
calfin_mat<-as.matrix(calfin)
sum(calfin_mat==calfin)
prod(dim(calfin))
calfin<-calfin_mat
rm(calfin_mat)
year_cf<-unname(calfin[1,2:(dim(calfin)[2])])
testthat::expect_equal(year,year_cf)
site_cf<-calfin[2:(dim(calfin)[1]),1]
testthat::expect_equal(site,site_cf)
calfin<-calfin[2:(dim(calfin)[1]),2:(dim(calfin)[2])]
testthat::expect_equal(dim(calfin),dim(pci))
rm(year_cf,site_cf)

#load the growing season temperature data
temp<-read.csv(file=paste0(datloc,"temp3to9.csv"),header=FALSE)
temp_mat<-as.matrix(temp)
sum(temp_mat==temp)
prod(dim(temp))
temp<-temp_mat
rm(temp_mat)
year_t<-unname(temp[1,2:(dim(temp)[2])])
testthat::expect_equal(year,year_t)
site_t<-temp[2:(dim(temp)[1]),1]
testthat::expect_equal(site,site_t)
temp<-temp[2:(dim(temp)[1]),2:(dim(temp)[2])]
testthat::expect_equal(dim(temp),dim(pci))
rm(year_t,site_t)

#These data just loaded have not been Box-Cox'd, or variance standardized, or detrended. They 
#are just the "raw" annual values, with gaps filled as described in the documenation for the 
#data and in the Plos Comp Biol paper.

#At this point, the level of cleaning (i.e., none so far except what was done previously
#for Plos Comp Biol) is similar to what was accomplished for kelp by 
#Kelp_AnalysisStep_CleanData_SuperBasic. Next step is to look at the next level
#of cleaning for kelp and do something appropriately similar for plankton



