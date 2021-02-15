#For computations associated with some spectral modelling math I did beginning on 2020 02 01.
#This generalizes (through functionization and other steps) the codes in 
#Approach_ModellingSpectral_MonterreyToMorro.R and Approach_ModellingSpectral_SantaBarbaraArea.R

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): parallel, wsyn, ncf
source("FillSeasonalMedian.R")
source("FindZeroNAYears.R")
source("SpectralTools.R")
source("AvgOffDiags.R")

#***
#Location for storing results and other prep
#***

resloc<-"../Results/Approach_ModellingSpectral_General/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

theseed<-101
BiasVariance<-0.5

saveRDS(theseed,file=paste0(resloc,"theseed.Rds"))
saveRDS(BiasVariance,file=paste0(resloc,"BiasVariance.Rds"))

#***
#Load the quarterly data with only super basic cleaning done to it
#***

#load the data
datloc<-"../Results/DataAfterSuperBasicCleaning/"
kelpQ<-readRDS(paste0(datloc,"Kelp_Quarterly_CleanedSuperBasic.Rds")) 
NO3Q<-readRDS(paste0(datloc,"NO3_Quarterly_CleanedSuperBasic.Rds")) 
wavesQ<-readRDS(paste0(datloc,"wave_Quarterly_CleanedSuperBasic.Rds")) 
locs<-readRDS(paste0(datloc,"Locs_CleanedSuperBasic.Rds"))
quarters<-readRDS(paste0(datloc,"Quarters_CleanedSuperBasic.Rds"))
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
#More extensive data prep - detrend, demean and variance standardize all time series individually
#***

kelp<-wsyn::cleandat(kelpQ,times=timesQ,clev=3)$cdat
NO3<-wsyn::cleandat(NO3Q,times=timesQ,clev=3)$cdat
waves<-wsyn::cleandat(wavesQ,times=timesQ,clev=3)$cdat

times<-timesQ

rm(kelpQ,wavesQ,NO3Q)
  
#***
#Which site is closest to pt conception, same for SB and some other locations, for later use
#***

#from google maps
MontBLatLon<-c(36.62393626182818, -121.8511283749223)
CarBLatLon<-c(36.52316357182972, -121.95535513873469)
MorBLatLon<-c(35.427553615598676, -120.89133343717805)
PtConcLatLon<-c(34.44830424836278,-120.47125509005721)
SBLatLon<-c(34.407876790966135, -119.68428533917637)
OxLatLon<-c(34.2310543912163, -119.2726684255394)
LALatLon<-c(34.00907538429378, -118.50422514984244)

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
#The main function for doing the analysis on a particular of locations
#***

#Uses the spectral mathematics from my notes begun on 2021 02 01 to get information about fractions of 
#synchrony explained by NO3, waves, and interactions. Assumes all the variables created above are in the 
#workspace when run.
#
#Args - a list with these named elements
#locstouse        Indices of locations to use. These refer to rows of the kelp, NO3 and waves data matrices above
#lags             Max regression lags to use, kelp then NO3 then waves. In quarters. So if you want the model to
#                   allow a kelp lag effect of 4 quarters and NO3 and wave lag effects from the current quarter 
#                   and last quarter, you input lags=c(4,1,1). Lags used will be 1 up to the value input for kelp,
#                   and 0 up to the value input for NO3, and the same for waves.
#frg              The frequency range to use for display
#fnpre            A prefix for file names of pdf plots exported by the function, with path but without the
#                   extension "pdf"
#
do_analysis<-function(Args)
{
  #**argument unpacking and error catching - global vars accessed here (kelp)
  if (class(Args)!="list" || length(Args)!=4)
  {
    stop("Error in do_analysis: Args must be a list of length 4")
  }
  if (any(names(Args)!=c("locstouse","lags","frg","fnpre")))
  {
    stop("Error in do_analysis: bad value for Args")
  }
  locstouse<-Args$locstouse
  lags<-Args$lags
  frg<-Args$frg
  fnpre<-Args$fnpre
  if (any(locstouse<1) || any(locstouse>dim(kelp)[1]))
  {
    stop("Error in do_analysis: bad value for lags")
  }
  if (length(locstouse)<length(unique(locstouse)))
  {
    stop("Error in do_analysis: repeat values not allowed in locstouse")
  }
  if (length(lags)!=3 ||
      any(lags<0) ||
      lags[1]<1)
  {
    stop("Error in do_analysis: bad value for lags")
  }
  
  #**filter the data to only the specified locations - global vars accessed here (kelp, NO3, waves)
  kelp_l<-kelp[locstouse,]
  NO3_l<-NO3[locstouse,]
  waves_l<-waves[locstouse,]
  numlocs_l<-length(locstouse)
  
  #**do regression to get coefficients for fB, fP1, fP2
  
  #construct the regression formula object
  form_rhs_terms<-c(paste0("kelp_l",1:lags[1]),paste0("NO3_l",0:lags[2]),paste0("waves_l",0:lags[2]))
  form<-as.formula(paste0("kelp_l0~",paste0(form_rhs_terms,collapse="+"),sep=""))
  
  #construct the regression data frame
  maxlag<-max(lags)
  regdat<-as.data.frame(matrix(NA,(dim(kelp_l)[1])*(dim(kelp_l)[2]-maxlag),length(form_rhs_terms)+1))
  names(regdat)<-c("kelp_l0",form_rhs_terms)
  for (counter in 0:lags[1]) #counter is the kelp lag here
  {
    regdat[,counter+1]<-as.vector(kelp_l[,(maxlag+1-counter):(lents-counter)])
  }
  for (counter in 0:lags[2]) #counter is the NO3 lag here
  {
    regdat[,counter+lags[1]+2]<-as.vector(NO3_l[,(maxlag+1-counter):(lents-counter)])
  }
  for (counter in 0:lags[3]) #counter is the waves lag here
  {
    regdat[,counter+lags[1]+1+lags[2]+1+1]<-as.vector(waves_l[,(maxlag+1-counter):(lents-counter)])
  }
  
  #do the regression, extract coefficients
  mod<-lm(formula=form,data=regdat)
  cm<-coef(mod)
  
  #**get residuals, cut data to common length, and compute the necessary spectral matrices
  
  #cut down the data - global vars accessed here (lents)
  delta<-matrix(unname(residuals(mod)),numlocs_l,lents-maxlag)
  kelp_l<-kelp_l[,(maxlag+1):lents]
  NO3_l<-NO3_l[,(maxlag+1):lents]
  waves_l<-waves_l[,(maxlag+1):lents]
  
  #get spectral matrices - global vars accessed here (BiasVariance)
  #e1 is NO3, e2 is waves, e is combined NO3 and then waves, w is kelp
  Sww<-myspecmatbrill(kelp_l,detrend=FALSE,cutsym=TRUE,forvar=FALSE,BiasVariance=BiasVariance)
  freq<-Sww$freq
  Sww<-Sww$spec
  S<-myspecmatbrill(rbind(NO3_l,waves_l,delta),detrend=FALSE,cutsym=TRUE,forvar=FALSE,BiasVariance=BiasVariance)
  S<-S$spec
  Se1e1<-S[1:numlocs_l,1:numlocs_l,]
  Se2e2<-S[(numlocs_l+1):(2*numlocs_l),(numlocs_l+1):(2*numlocs_l),]
  Se1e2<-S[1:numlocs_l,(numlocs_l+1):(2*numlocs_l),]
  Se2e1<-S[(numlocs_l+1):(2*numlocs_l),1:numlocs_l,]
  Sdd<-S[(2*numlocs_l+1):(3*numlocs_l),(2*numlocs_l+1):(3*numlocs_l),]
  Se1d<-S[1:numlocs_l,(2*numlocs_l+1):(3*numlocs_l),]
  Sde1<-S[(2*numlocs_l+1):(3*numlocs_l),1:numlocs_l,]
  Se2d<-S[(numlocs_l+1):(2*numlocs_l),(2*numlocs_l+1):(3*numlocs_l),]
  Sde2<-S[(2*numlocs_l+1):(3*numlocs_l),(numlocs_l+1):(2*numlocs_l),]
  
  #**construct fB, fP1, fP2
  mu<-exp(-2*pi*complex(real=0,imaginary=1)*freq)
  fP1<-complex(real=numeric(length(freq)),imaginary=numeric(length(freq)))
  for (counter in 0:lags[2])
  {
    fP1<-fP1+unname(cm[lags[1]+2+counter])*mu^counter
  }
  fP2<-complex(real=numeric(length(freq)),imaginary=numeric(length(freq)))
  for (counter in 0:lags[3])
  {
    fP2<-fP2+unname(cm[lags[1]+1+lags[2]+1+1+counter])*mu^counter
  }
  fB<-complex(real=numeric(length(freq))+1,imaginary=numeric(length(freq)))
  for (counter in 1:lags[1])
  {
    fB<-fB-unname(cm[counter+1])*mu^counter
  }
  
  #**get the terms of the spectral equation (p. 6 of my math notes)
  
  #term 1
  T1pre<-((Mod(fP1))^2)/((Mod(fB))^2) 
  T1<-aperm(array(T1pre,c(dim(Se1e1)[3],dim(Se1e1)[1:2])),c(2,3,1))
  T1<-T1*Se1e1
  
  #term 2
  T2pre<-((Mod(fP2))^2)/((Mod(fB))^2) 
  T2<-aperm(array(T2pre,c(dim(Se2e2)[3],dim(Se2e2)[1:2])),c(2,3,1))
  T2<-T2*Se2e2
  
  #term 3
  T3pre<-1/((Mod(fB))^2) 
  T3<-aperm(array(T3pre,c(dim(Sdd)[3],dim(Sdd)[1:2])),c(2,3,1))
  T3<-T3*Sdd
  
  #term 4
  T4pre1<-fP1*Conj(fP2)/((Mod(fB))^2)
  T4pre2<-fP2*Conj(fP1)/((Mod(fB))^2)
  T4p1<-aperm(array(T4pre1,c(dim(Se1e2)[3],dim(Se1e2)[1:2])),c(2,3,1))
  T4p1<-T4p1*Se1e2
  T4p2<-aperm(array(T4pre2,c(dim(Se2e1)[3],dim(Se2e1)[1:2])),c(2,3,1))
  T4p2<-T4p2*Se2e1
  T4<-T4p1+T4p2  

  #term 5
  T5pre1<-fP1/((Mod(fB))^2)
  T5pre2<-Conj(fP1)/((Mod(fB))^2)
  T5p1<-aperm(array(T5pre1,c(dim(Se1d)[3],dim(Se1d)[1:2])),c(2,3,1))
  T5p1<-T5p1*Se1d
  T5p2<-aperm(array(T5pre2,c(dim(Sde1)[3],dim(Sde1)[1:2])),c(2,3,1))
  T5p2<-T5p2*Sde1
  T5<-T5p1+T5p2
  
  #term 6
  T6pre1<-fP2/((Mod(fB))^2)
  T6pre2<-Conj(fP2)/((Mod(fB))^2)
  T6p1<-aperm(array(T6pre1,c(dim(Se2d)[3],dim(Se2d)[1:2])),c(2,3,1))
  T6p1<-T6p1*Se2d
  T6p2<-aperm(array(T6pre2,c(dim(Sde2)[3],dim(Sde2)[1:2])),c(2,3,1))
  T6p2<-T6p2*Sde2
  T6<-T6p1+T6p2
  
  #**average real parts across pairs of distinct locations
  totsync<-avgoffdiags(Re(Sww))
  T1_avg<-avgoffdiags(Re(T1))
  T2_avg<-avgoffdiags(Re(T2))
  T3_avg<-avgoffdiags(Re(T3))
  T4_avg<-avgoffdiags(Re(T4))
  T5_avg<-avgoffdiags(Re(T5))
  T6_avg<-avgoffdiags(Re(T6))
  
  #**get various statistics to return
  
  #linear model stuff
  res<-c()
  res[1]<-summary(mod)$r.squared
  names(res)[1]<-"model r sq"
  res[2:(length(cm))]<-cm[2:length(cm)]
  names(res)[2:(length(cm))]<-paste0("coef_",names(cm)[2:(length(cm))])
  
  #annual timescale stuff
  h<-c()
  inds<-which(freq>.2 & freq<=.3)
  h[1]<-mean(T1_avg[inds])  
  h[2]<-mean(T2_avg[inds])
  h[3]<-mean(T4_avg[inds])
  h[4]<-mean(totsync[inds])
  names(h)<-c("Ann sync NO3","Ann sync waves","Ann sync NO3/waves","Ann sync tot")
  res<-c(res,h)
    
  #2-4 years stuff
  h<-c()
  inds<-which(freq>1/16 & freq<=1/8)
  h[1]<-mean(T1_avg[inds])  
  h[2]<-mean(T2_avg[inds])
  h[3]<-mean(T4_avg[inds])
  h[4]<-mean(totsync[inds])
  names(h)<-c("2to4yr sync NO3","2to4yr sync waves","2to4yr sync NO3/waves","2to4yr sync tot")
  res<-c(res,h)
  
  #>4 years stuff
  h<-c()
  inds<-which(freq<=1/16)
  h[1]<-mean(T1_avg[inds])  
  h[2]<-mean(T2_avg[inds])
  h[3]<-mean(T4_avg[inds])
  h[4]<-mean(totsync[inds])
  names(h)<-c(">4yr sync NO3",">4yr sync waves",">4yr sync NO3/waves",">4yr sync tot")
  res<-c(res,h)

  #basics
  h<-c(min(locstouse),max(locstouse),numlocs_l)
  names(h)<-c("first loc","last loc","num locs")
  res<-c(res,h)
  
  #**constraint to certain frequencies, for plotting
  totsync<-totsync[frg[1]<freq & freq<=frg[2]]
  T1_avg<-T1_avg[frg[1]<freq & freq<=frg[2]]
  T2_avg<-T2_avg[frg[1]<freq & freq<=frg[2]]
  T3_avg<-T3_avg[frg[1]<freq & freq<=frg[2]]
  T4_avg<-T4_avg[frg[1]<freq & freq<=frg[2]]
  T5_avg<-T5_avg[frg[1]<freq & freq<=frg[2]]
  T6_avg<-T6_avg[frg[1]<freq & freq<=frg[2]]
  freq<-freq[frg[1]<freq & freq<=frg[2]]
  
  #**make and export the plot
  jpeg(paste0(resloc,fnpre,"_MainPlot.jpg"),quality=95)
  sumexpl<-T1_avg+T2_avg+T4_avg
  sumunexpl<-T3_avg+T5_avg+T6_avg
  sumterms<-sumexpl+sumunexpl
  ylimits<-range(totsync,T1_avg,T2_avg,T4_avg,sumexpl,sumterms)
  plot(freq,totsync,type="l",xlab="Frequency, cycles per quarter",ylab="Synchrony contribution",col="black",ylim=ylimits)
  lines(freq,T1_avg,type="l",col="green") #direct effects of NO3
  lines(freq,T2_avg,type="l",col="blue") #direct effects of waves
  lines(freq,T4_avg,type="l",lty="solid",col="green") 
  lines(freq,T4_avg,type="l",lty="dashed",col="blue") #interactions between waves and NO3
  lines(freq,sumexpl,type="l",col="yellow") #explained stuff - sum of all the terms relating only to NO3 and waves (direct NO3 and wave effects, and interactions between those two)
  lines(freq,sumterms,type="l",lty="dashed") #total of all terms, should approx equal the plot of totsync
  lines(freq,sumunexpl,type="l",lty="solid",col="red") #unexplained
  lines(freq,rep(0,length(freq)))
  lines(rep(1/4,2),ylimits) #plot a vertical line at the annual frequency
  lines(rep(1/8,2),ylimits) #plot a vertical line at the biannual frequency
  lines(rep(1/16,2),ylimits) #plot another one at the 4-year timescale
  dev.off()
  
  #return results
  return(res)
}

#***
#some tests of the function, based on comparison to analyses done elsewhere using seperate code
#***

Args<-list(locstouse=25:220,lags=c(4,1,1),frg=c(0,.3),fnpre="Test1")
do_analysis(Args)
Args<-list(locstouse=25:220,lags=c(8,1,1),frg=c(0,.3),fnpre="Test2")
do_analysis(Args)

#***
#Do a big run of the function for sets of sites up and down the coast
#***

#assemble a list of argument lists for a call to mclapply, and do that call
alldists<-ncf::gcdist(locs$Lon,locs$Lat)
allargs<-list()
allargsind<-1
centerinds<-c()
for (counter in 1:(dim(locs)[1]))
{
  locstouse<-which(alldists[counter,]<25)
  if (length(locstouse)<20){ next }
  centerinds<-c(centerinds,counter)
  curargs1<-list(locstouse=locstouse,lags=c(4,1,1),frg=c(0,.3),fnpre=paste0("KelpLag_4_Center_",counter))
  curargs2<-list(locstouse=locstouse,lags=c(8,1,1),frg=c(0,.3),fnpre=paste0("KelpLag_8_Center_",counter))
  allargs[[allargsind]]<-curargs1
  allargs[[allargsind+1]]<-curargs2
  allargsind<-allargsind+2
}
allres<-parallel::mclapply(FUN=do_analysis,X=allargs,mc.cores=10)

saveRDS(allargs,file=paste0(resloc,"allargs.Rds"))
saveRDS(allres,file=paste0(resloc,"allres.Rds"))

#***
#now organize summary results 
#***

#put the results in a data frame for lag 4
allres4<-matrix(NA,length(allres)/2,length(allres[[1]]))
colnames(allres4)<-names(allres[[1]])
for (counter in seq(from=1,by=2,to=length(allres)))
{
  allres4[(counter+1)/2,]<-allres[[counter]]
}
saveRDS(allres4,paste0(resloc,"allres4.Rds"))

#same for lag 8
allres8<-matrix(NA,length(allres)/2,length(allres[[2]]))
colnames(allres8)<-names(allres[[2]])
for (counter in seq(from=1,by=2,to=length(allres)))
{
  allres8[(counter+1)/2,]<-allres[[counter+1]]
}
saveRDS(allres8,paste0(resloc,"allres8.Rds"))

#get variables we can use as x axes of the plots below
mdlocind<-c() #will be median of indices of locations used in each run
mnloclat<-c() #will be mean of latitudes of locations used
mnloclon<-c() #will be mean of longitudes of locations used
for (counter in seq(from=1,by=2,to=length(allargs)))
{
  mdlocind[(counter+1)/2]<-median(allargs[[counter]]$locstouse)
  mnloclat[(counter+1)/2]<-mean(locs$Lat[allargs[[counter]]$locstouse])
  mnloclon[(counter+1)/2]<-mean(locs$Lon[allargs[[counter]]$locstouse])
}
ctrloclat<-locs$Lat[centerinds]
ctrloclon<-locs$Lon[centerinds]

#***
#make plots which can (hopefully) be interpreted to get a sense for whether results are more or less
#consistent across central CA, and more or less consistent across southern CA
#***

put_loc_lines<-function()
{
  text(PtConcInd,ylimits[2],"Pt. Conc.",srt=90,adj=c(1,0))
  lines(rep(PtConcInd,2),ylimits)
  
  text(MontBInd,ylimits[2],"Mont. Bay",srt=90,adj=c(1,0))
  lines(rep(MontBInd,2),ylimits,lty="dashed")
  
  text(CarBInd,ylimits[2],"Car. Bay",srt=90,adj=c(1,0))
  lines(rep(CarBInd,2),ylimits,lty="dashed")
  
  text(MorBInd,ylimits[2],"Mor. Bay",srt=90,adj=c(1,0))
  lines(rep(MorBInd,2),ylimits,lty="dashed")
  
  text(SBInd,ylimits[2],"SB",srt=90,adj=c(1,0))
  lines(rep(SBInd,2),ylimits,lty="dashed")
    
  lines(rep(OxInd,2),ylimits,lty="dashed")
  text(OxInd,ylimits[2],"Oxnard",srt=90,adj=c(1,0))
  
  lines(rep(LAInd,2),ylimits,lty="dashed")
  text(LAInd,ylimits[2],"LA",srt=90,adj=c(1,0))
}

#make plots about R^2 - how well does the ARMA capture the dynamics
pdf(paste0(resloc,"LinearModel_Rsq_v_centerind.pdf"))
ylimits<-range(allres4[,'model r sq'])
plot(centerinds,allres4[,'model r sq'],type='b',xlab="Center index",ylab="Linear model R sq")
put_loc_lines()
dev.off()

#make plots about coefficients - kelp coefficients first
pdf(paste0(resloc,"LinearModel_KelpCoefs_v_centerind.pdf"))
ylimits<-range(allres4[,'coef_kelp_l1'],allres4[,'coef_kelp_l2'],allres4[,'coef_kelp_l2'],allres4[,'coef_kelp_l4'])
plot(centerinds,allres4[,'coef_kelp_l1'],type='b',xlab="Center index",ylab="Kelp coefficient",pch=1,
     ylim=ylimits)
lines(centerinds,allres4[,'coef_kelp_l2'],type='b',pch=2)
lines(centerinds,allres4[,'coef_kelp_l3'],type='b',pch=22)
lines(centerinds,allres4[,'coef_kelp_l4'],type='b',pch=11)
put_loc_lines()
dev.off()

#now NO3 coefficients
pdf(paste0(resloc,"LinearModel_NO3Coefs_v_centerind.pdf"))
ylimits<-range(allres4[,'coef_NO3_l0'],allres4[,'coef_NO3_l1'])
plot(centerinds,allres4[,'coef_NO3_l0'],type='b',xlab="Center index",ylab="NO3 coefficient",pch=1,
     ylim=ylimits)
lines(centerinds,allres4[,'coef_NO3_l1'],type='b',pch=2)
put_loc_lines()
dev.off()

#now waves coefficients
pdf(paste0(resloc,"LinearModel_wavesCoefs_v_centerind.pdf"))
ylimits<-range(allres4[,'coef_waves_l0'],allres4[,'coef_waves_l1'])
plot(centerinds,allres4[,'coef_waves_l0'],type='b',xlab="Center index",ylab="Waves coefficient",pch=1,
     ylim=ylimits)
lines(centerinds,allres4[,'coef_waves_l1'],type='b',pch=2)
put_loc_lines()
dev.off()

#now plot info about annual-timescale sycnhrony
pdf(paste0(resloc,"Synchrony_AnnualTimescale_v_centerind.pdf"))
ylimits<-range(0,allres4[,'Ann sync tot'],allres4[,'Ann sync NO3'],allres4[,'Ann sync waves'],allres4[,'Ann sync NO3/waves'])
plot(centerinds,allres4[,'Ann sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, annual timescales",pch=20,
     ylim=ylimits)
lines(centerinds,allres4[,'Ann sync NO3'],type='b',pch=20,col="green")
lines(centerinds,allres4[,'Ann sync waves'],type='b',pch=20,col="blue")
lines(centerinds,allres4[,'Ann sync NO3/waves'],type='b',pch=19,col="green")
lines(centerinds,allres4[,'Ann sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
lines(centerinds,allres4[,'Ann sync NO3']+allres4[,'Ann sync waves']+allres4[,'Ann sync NO3/waves'],
      type="b",col="yellow",pch=20)
lines(centerinds,rep(0,length(centerinds)),type="l")
put_loc_lines()
dev.off()

#***DAN: IDEA, maybe the reason why synchrony is so much less in SoCal is because waves are much
#less synchronous, *and* their effects are no longer aligned with the effects of NO3. That means
#you lose both the wave and NO3/wave interaction effects on synch. Note the diff between the black and
#yellow lines in the above plot appears pretty similar across all of CA, indicating we are explaining
#a similar fraction of synchrony. This is consistent with the hypothesis.

#to follow up on the above idea, plot the difference
pdf(paste0(resloc,"SynchronyDiff_AnnualTimescale_TotMinusExplained_v_centerind.pdf"))
plot(centerinds,allres4[,'Ann sync tot']-(allres4[,'Ann sync NO3']+allres4[,'Ann sync waves']+allres4[,'Ann sync NO3/waves']),type='b')
put_loc_lines()
dev.off()
h<-allres4[,'Ann sync tot']-(allres4[,'Ann sync NO3']+allres4[,'Ann sync waves']+allres4[,'Ann sync NO3/waves'])
mean(h[CarBInd:PtConcInd])
mean(h[PtConcInd:SBInd]) 
#these are not wildly different, which may support the hypothesis.

#now do 2-4 year timescales
pdf(paste0(resloc,"Synchrony_2to4yrTimescales_v_centerind.pdf"))
ylimits<-range(0,allres4[,'2to4yr sync tot'],allres4[,'2to4yr sync NO3'],allres4[,'2to4yr sync waves'],allres4[,'2to4yr sync NO3/waves'])
plot(centerinds,allres4[,'2to4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, 2-4 yr timescales",pch=20,
     ylim=ylimits)
lines(centerinds,allres4[,'2to4yr sync NO3'],type='b',pch=20,col="green")
lines(centerinds,allres4[,'2to4yr sync waves'],type='b',pch=20,col="blue")
lines(centerinds,allres4[,'2to4yr sync NO3/waves'],type='b',pch=19,col="green")
lines(centerinds,allres4[,'2to4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
lines(centerinds,allres4[,'2to4yr sync NO3']+allres4[,'2to4yr sync waves']+allres4[,'2to4yr sync NO3/waves'],
      type="b",col="yellow",pch=20)
lines(centerinds,rep(0,length(centerinds)),type="l")
put_loc_lines()
dev.off()

#now do >4 yr timescales
pdf(paste0(resloc,"Synchrony_gt4yrTimescales_v_centerind.pdf"))
ylimits<-range(0,allres4[,'>4yr sync tot'],allres4[,'>4yr sync NO3'],allres4[,'>4yr sync waves'],allres4[,'>4yr sync NO3/waves'])
plot(centerinds,allres4[,'>4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, >4 yr timescales",pch=20,
     ylim=ylimits)
lines(centerinds,allres4[,'>4yr sync NO3'],type='b',pch=20,col="green")
lines(centerinds,allres4[,'>4yr sync waves'],type='b',pch=20,col="blue")
lines(centerinds,allres4[,'>4yr sync NO3/waves'],type='b',pch=19,col="green")
lines(centerinds,allres4[,'>4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
lines(centerinds,allres4[,'>4yr sync NO3']+allres4[,'>4yr sync waves']+allres4[,'>4yr sync NO3/waves'],
      type="b",col="yellow",pch=20)
lines(centerinds,rep(0,length(centerinds)),type="l")
put_loc_lines()
dev.off()

#plot numbers of sites used 
pdf(paste0(resloc,"Numsitesused_v_centerind.pdf"))
ylimits<-range(allres4[,'num locs'])
plot(centerinds,allres4[,'num locs'],type='b',xlab="Center index",ylab="Number of sites used")
lines(rep(PtConcInd,2),ylimits)
lines(rep(SBInd,2),ylimits,lty="dashed")
put_loc_lines()
dev.off()

#***DAN: THOUGHTS: I think if I choose a region in Central CA and one in SoCal (centered slightly west 
#of SB) and contrast them, that will be good. I think generally results are roughly consistent across 
#Central CA and across SoCal from Pt Conception to Oxnard (and after that the sites are too spaced out 
#to do this analsis). The idea would be to present one analysis from Central Cal and one from So Cal 
#and then say those are representative of if we had used different areas. Maybe do two from Central CA, 
#since it's bigger than the usable So Cal area. Then you would not have to say your selected area is
#representative, since you'd be using the bulk of the Central CA sites. So split up Carmel Bay to Morro 
#Bay into two sets of sites, and use Pt Conc to Oxnard for the So Cal sites.

#***
#Select some regions
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

#So my three roughly comparable comparable regions are SBlocstouse, CC1locstouse, CC2locstouse

#***
#Now do the analysis for these three regions and save
#***

Args<-list(locstouse=CC1locstouse,lags=c(4,1,1),frg=c(0,.3),fnpre="RegionalAnalysis_CentCal1_KelpLag4")
res_CentCal1_KelpLag4<-do_analysis(Args)
saveRDS(res_CentCal1_KelpLag4,paste0(resloc,"res_CentCal1_KelpLag4.Rds"))

Args<-list(locstouse=CC1locstouse,lags=c(8,1,1),frg=c(0,.3),fnpre="RegionalAnalysis_CentCal1_KelpLag8")
res_CentCal1_KelpLag8<-do_analysis(Args)
saveRDS(res_CentCal1_KelpLag8,paste0(resloc,"res_CentCal1_KelpLag8.Rds"))

Args<-list(locstouse=CC2locstouse,lags=c(4,1,1),frg=c(0,.3),fnpre="RegionalAnalysis_CentCal2_KelpLag4")
res_CentCal2_KelpLag4<-do_analysis(Args)
saveRDS(res_CentCal2_KelpLag4,paste0(resloc,"res_CentCal2_KelpLag4.Rds"))

Args<-list(locstouse=CC2locstouse,lags=c(8,1,1),frg=c(0,.3),fnpre="RegionalAnalysis_CentCal2_KelpLag8")
res_CentCal2_KelpLag8<-do_analysis(Args)
saveRDS(res_CentCal2_KelpLag8,paste0(resloc,"res_CentCal2_KelpLag8.Rds"))

Args<-list(locstouse=SBlocstouse,lags=c(4,1,1),frg=c(0,.3),fnpre="RegionalAnalysis_SoCal_KelpLag4")
res_SoCal_KelpLag4<-do_analysis(Args)
saveRDS(res_SoCal_KelpLag4,paste0(resloc,"res_SoCal_KelpLag4.Rds"))

Args<-list(locstouse=SBlocstouse,lags=c(8,1,1),frg=c(0,.3),fnpre="RegionalAnalysis_SoCal_KelpLag8")
res_SoCal_KelpLag8<-do_analysis(Args)
saveRDS(res_SoCal_KelpLag8,paste0(resloc,"res_SoCal_KelpLag8.Rds"))

#***DAN: IDEA: Change do_analysis so it takes an argument which controls the format of the plot generated
#(jpg or pdf), so that you can make the above plots pdf while still making the many plots up and down the
#coast jpgs, for easy rapid viewing.
#IDEA: Find another way to plot that does not excessively de-emphasize the >4yr timescale range. Maybe
#plot against timescale, though that may de-emphasize the annual timescale. Maybe do line is done
#with wavelets, which uses log-spaced plotting - review it and consider making the change here.
#IDEA: Consider plotting the y-axis differently so that the huge annual peak does not drown out what is
#happening at >4 yr timescales. Maybe plot two panels, one for annual, one for super-annual timescales,
#with different y axis extents?