#For computations associated with some spectral modelling math I did beginning on 2020 02 01.
#This generalizes (through functionization and other steps) the codes in 
#Approach_ModellingSpectral_MonterreyToMorro.R and Approach_ModellingSpectral_SantaBarbaraArea.R

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): parallel, wsyn, ncf, graphics, latex2exp
source("SpectralTools.R")
source("AvgOffDiags.R")
source("PlotPhaseFunc.R")
source("WhichBreak.R")
source("plotColLine.R")

#***
#Location for storing results and other prep
#***

resloc<-"../Results/Kelp_SpectralDecompResults/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

theseed<-101
BiasVariance<-0.5

saveRDS(theseed,file=paste0(resloc,"theseed.Rds"))
saveRDS(BiasVariance,file=paste0(resloc,"BiasVariance.Rds"))

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
#                   extension "pdf" of "jpg"
#plottype         Use "jpg" to make jpg plots, "pdf" to make pdf plots
#PlotForMSArgs    Formatting arguments used for generating final, presentation-quality plots. If this is missing it 
#                   means don't generate those plots.
#
do_analysis<-function(Args)
{
  #**argument unpacking and error catching - global vars accessed here (kelp)
  if (class(Args)!="list" || !(length(Args)%in%c(5,6)))
  {
    stop("Error in do_analysis: Args must be a list of length 5 or 6")
  }
  if (any(names(Args)[1:5]!=c("locstouse","lags","frg","fnpre","plottype")))
  {
    stop("Error in do_analysis: bad value for Args")
  }
  locstouse<-Args$locstouse
  lags<-Args$lags
  frg<-Args$frg
  fnpre<-Args$fnpre
  plottype<-Args$plottype
  if (length(Args)==6)
  {
    PlotForMSArgs<-Args[[6]]
  }else
  {
    PlotForMSArgs<-NA
  }
  if (any(locstouse<1) || any(locstouse>dim(kelp)[1]))
  {
    stop("Error in do_analysis: bad value for locstouse")
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
  if (!(plottype %in% c("jpg","pdf")))
  {
    stop("Error in do_analysis: bad value for plottype")
  }
  
  #**filter the data to only the specified locations - global vars accessed here (kelp, NO3, waves)
  kelp_l<-kelp[locstouse,]
  NO3_l<-NO3[locstouse,]
  waves_l<-waves[locstouse,]
  numlocs_l<-length(locstouse)
  
  #**get variances of kelp times series, which we'll need later and want to base on the whole kelp
  #time series instead of on the truncated version which we're about to restrict to for other purposes
  kelpvars<-apply(FUN=var,X=kelp_l,MARGIN=1)
  
  #**do regression to get coefficients for fB, fP1, fP2
  
  #construct the regression formula object
  form_rhs_terms<-c(paste0("kelp_l",1:lags[1]),paste0("NO3_l",0:lags[2]),paste0("waves_l",0:lags[3]))
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
  
  #**average real parts across pairs of distinct locations, WITHOUT normalizing by kelp variances
  totsync_cov<-avgoffdiags(Re(Sww))
  T1_avg_cov<-avgoffdiags(Re(T1))
  T2_avg_cov<-avgoffdiags(Re(T2))
  T3_avg_cov<-avgoffdiags(Re(T3))
  T4_avg_cov<-avgoffdiags(Re(T4))
  T5_avg_cov<-avgoffdiags(Re(T5))
  T6_avg_cov<-avgoffdiags(Re(T6))
  Se1e1_avg_cov<-avgoffdiags(Re(Se1e1))
  Se2e2_avg_cov<-avgoffdiags(Re(Se2e2))
  
  Se1e2_avg_cov<-avgoffdiags(Se1e2)
  
  #**same thing WITH normalizing by kelp variances
  totsync_cor<-avgoffdiags(Re(Sww),kelpvars)
  T1_avg_cor<-avgoffdiags(Re(T1),kelpvars)
  T2_avg_cor<-avgoffdiags(Re(T2),kelpvars)
  T3_avg_cor<-avgoffdiags(Re(T3),kelpvars)
  T4_avg_cor<-avgoffdiags(Re(T4),kelpvars)
  T5_avg_cor<-avgoffdiags(Re(T5),kelpvars)
  T6_avg_cor<-avgoffdiags(Re(T6),kelpvars)
  Se1e1_avg_cor<-avgoffdiags(Re(Se1e1),kelpvars)
  Se2e2_avg_cor<-avgoffdiags(Re(Se2e2),kelpvars)
  
  #**get various statistics to return
  
  #linear model stuff
  res<-c()
  res[1]<-summary(mod)$r.squared
  names(res)[1]<-"model r sq"
  res[2:(length(cm))]<-cm[2:length(cm)]
  names(res)[2:(length(cm))]<-paste0("coef_",names(cm)[2:(length(cm))])

  #timescale-specific synchrony quantities, gathered by a utility function, "get_stats", which is below
  stats_cov<-get_stats(freq,T1_avg_cov,T2_avg_cov,T4_avg_cov,totsync_cov)
  stats_cor<-get_stats(freq,T1_avg_cor,T2_avg_cor,T4_avg_cor,totsync_cor)
  names(stats_cov)<-paste0("cov_",names(stats_cov))
  names(stats_cor)<-paste0("cor_",names(stats_cor))
  res<-c(res,stats_cov,stats_cor)

  #basics stats
  h<-c(min(locstouse),max(locstouse),numlocs_l)
  names(h)<-c("first loc","last loc","num locs")
  res<-c(res,h)
  
  #**make plots - this is the "_MainPlot" plot that shows all the components of synchrony
  make_plots(freq,frg,totsync_cov,T1_avg_cov,T2_avg_cov,T3_avg_cov,T4_avg_cov,T5_avg_cov,T6_avg_cov,resloc,
             fnpre=paste0(fnpre,"_cov"),plottype)
  make_plots(freq,frg,totsync_cor,T1_avg_cor,T2_avg_cor,T3_avg_cor,T4_avg_cor,T5_avg_cor,T6_avg_cor,resloc,
             fnpre=paste0(fnpre,"_cor"),plottype)
  
  #make plots of components fP1, fP2, etc., for understanding, with argument and phase on separate axes
  make_plot_component(freq=freq,frg=frg,comp=fP1,plotname=paste0(fnpre,"_fP1"),plottype=plottype)
  make_plot_component(freq=freq,frg=frg,comp=fP2,plotname=paste0(fnpre,"_fP2"),plottype=plottype)
  make_plot_component(freq=freq,frg=frg,comp=fP1*Conj(fP2),plotname=paste0(fnpre,"_fP1ConjfP2"),
                      plottype=plottype)
  
  #make plots of Se1e2 summed across all i \neq j
  h<-Se1e2
  for (counter in 1:(dim(h)[1]))
  {
    h[counter,counter,]<-NA
  }
  h<-apply(FUN=sum,X=h,MARGIN=3,na.rm=TRUE)
  make_plot_component(freq=freq,frg=frg,comp=h,plotname=paste0(fnpre,"_cov_Cross"),plottype=plottype)
  
  #plot this same thing times fP1*Conj(fP2)
  make_plot_component(freq=freq,frg=frg,comp=fP1*Conj(fP2)*h,plotname=paste0(fnpre,"_cov_All"),plottype=plottype)
  
  #now do the same things for the "cor" approach to synchrony
  h<-Se1e2
  for (counter in 1:(dim(h)[1]))
  {
    h[counter,counter,]<-NA
  }
  vivj<-outer(X=sqrt(kelpvars),Y=sqrt(kelpvars),FUN="*")
  for (counter in 1:(dim(h)[3]))
  {
    h[,,counter]<-h[,,counter]*vivj
  }
  h<-apply(FUN=sum,X=h,MARGIN=3,na.rm=TRUE)
  make_plot_component(freq=freq,frg=frg,comp=h,plotname=paste0(fnpre,"_cor_Cross"),plottype=plottype)
  
  #plot this same thing times fP1*Conj(fP2)
  make_plot_component(freq=freq,frg=frg,comp=fP1*Conj(fP2)*h,plotname=paste0(fnpre,"_cor_All"),plottype=plottype)
  
  #make plots of synchrony of NO3 and waves
  make_plots_noise_sync(freq=freq,frg=frg,S_avg=Se1e1_avg_cov,fnpre=paste0(fnpre,"_cov_NO3"),plottype=plottype)
  make_plots_noise_sync(freq=freq,frg=frg,S_avg=Se1e1_avg_cor,fnpre=paste0(fnpre,"_cor_NO3"),plottype=plottype)
  make_plots_noise_sync(freq=freq,frg=frg,S_avg=Se2e2_avg_cov,fnpre=paste0(fnpre,"_cov_Waves"),plottype=plottype)
  make_plots_noise_sync(freq=freq,frg=frg,S_avg=Se2e2_avg_cor,fnpre=paste0(fnpre,"_cor_Waves"),plottype=plottype)
  
  #make presentation-quality plots if desired
  if (class(PlotForMSArgs)=="list")
  {
    #make the plot that shows the decomposition of synchrony into main effects and interaction effects
    h<-list(MainPlot_PanLabs=PlotForMSArgs$cov_MainPlot_PanLabs,
            MainPlot_short_ylim=PlotForMSArgs$cov_MainPlot_short_ylim,
            MainPlot_mid_ylim=PlotForMSArgs$cov_MainPlot_mid_ylim,
            MainPlot_long_ylim=PlotForMSArgs$cov_MainPlot_long_ylim            )
    make_plots_ForMS(freq,frg,totsync_cov,T1_avg_cov,T2_avg_cov,T3_avg_cov,T4_avg_cov,T5_avg_cov,T6_avg_cov,resloc,
               fnpre=paste0(fnpre,"_cov"),plottype,PlotForMSArgs=h)

    #make the plot that explains the main effects
    explain_direct_Moran_ForMS(freq=freq,frg=frg,fP1=fP1,fP2=fP2,fB=fB,Se1e1_avg=Se1e1_avg_cov,Se2e2_avg=Se2e2_avg_cov,
                               fnpre=paste0(fnpre,"_cov"),plottype=plottype,PlotForMSArgs=PlotForMSArgs)
    
    #make the plot that explains the interactions between Moran effects
    explain_interacting_Moran(freq=freq,frg=frg,fP1=fP1,fP2=fP2,Se1e2_avg=Se1e2_avg_cov,
                              fnpre=paste0(fnpre,"_cov"),plottype=plottype,PlotForMSArgs=PlotForMSArgs)
  }
  
  return(res)
}

#A function for making the manuscript plots explaining interactions between Moran effects
#
explain_interacting_Moran<-function(freq,frg,fP1,fP2,Se1e2_avg,
                          fnpre,plottype,PlotForMSArgs)
{
  #unpack plotting params
  fP1_short_ylims<-PlotForMSArgs$ExplainInteractingMoran_fP1_short_ylims
  fP1_mid_ylims<-PlotForMSArgs$ExplainInteractingMoran_fP1_mid_ylims
  fP1_long_ylims<-PlotForMSArgs$ExplainInteractingMoran_fP1_long_ylims
  fP2_short_ylims<-PlotForMSArgs$ExplainInteractingMoran_fP2_short_ylims
  fP2_mid_ylims<-PlotForMSArgs$ExplainInteractingMoran_fP2_mid_ylims
  fP2_long_ylims<-PlotForMSArgs$ExplainInteractingMoran_fP2_long_ylims
  fP1ConjfP2_short_ylims<-PlotForMSArgs$ExplainInteractingMoran_fP1ConjfP2_short_ylims
  fP1ConjfP2_mid_ylims<-PlotForMSArgs$ExplainInteractingMoran_fP1ConjfP2_mid_ylims
  fP1ConjfP2_long_ylims<-PlotForMSArgs$ExplainInteractingMoran_fP1ConjfP2_long_ylims
  Se1e2_short_ylims<-PlotForMSArgs$ExplainInteractingMoran_Se1e2_short_ylims
  Se1e2_mid_ylims<-PlotForMSArgs$ExplainInteractingMoran_Se1e2_mid_ylims
  Se1e2_long_ylims<-PlotForMSArgs$ExplainInteractingMoran_Se1e2_long_ylims
  Prod_short_ylims<-PlotForMSArgs$ExplainInteractingMoran_Prod_short_ylims
  Prod_mid_ylims<-PlotForMSArgs$ExplainInteractingMoran_Prod_mid_ylims
  Prod_long_ylims<-PlotForMSArgs$ExplainInteractingMoran_Prod_long_ylims
  ExplainInteractingMoran_PanLabs<-PlotForMSArgs$ExplainInteractingMoran_PanLabs
  
  #constraint to certain frequencies, for plotting
  fP1<-fP1[frg[1]<freq & freq<=frg[2]]
  fP2<-fP2[frg[1]<freq & freq<=frg[2]]
  Se1e2_avg<-Se1e2_avg[frg[1]<freq & freq<=frg[2]]
  freq<-freq[frg[1]<freq & freq<=frg[2]]
  
  #we are going to plot against log2 timescale, on three panels for different timescale ranges, so
  #prep for that
  timescales<-1/freq
  l2timescales<-log2(timescales)
  xlimits<-range(l2timescales)
  xlimits_short<-c(xlimits[1],log2(8))
  xlimits_mid<-c(log2(8),log2(16))
  xlimits_long<-c(log2(16),xlimits[2])
  
  #panel layout quantities
  gap<-.15
  yaxnumwd<-.4
  yaxlabwd<-.4
  yaxwd<-yaxnumwd+yaxlabwd
  xaxht<-yaxwd
  panht<-1.25
  panwd_short<-2
  panwd_mid<-1
  panwd_long<-2
  totpanwd<-panwd_short+panwd_mid+panwd_long
  totht<-xaxht+5*panht+5*gap
  totwd<-yaxwd+2*yaxnumwd+totpanwd+gap
  lnwid<-2.5
  
  #get the plotting device
  if (plottype=="jpg")
  {
    jpeg(paste0(resloc,fnpre,"_ExplainInteractingMoran.jpg"),quality=95,width=totwd,height=totht,units="in",res=300)
  }
  if (plottype=="pdf")
  {
    pdf(paste0(resloc,fnpre,"_ExplainInteractingMoran.pdf"),width=totwd,height=totht)
  }
  
  #technical stuff for the circular colorbar for phases
  breaks<-seq(from=-pi,to=pi,length.out=201)
  colbarQuad1<-rgb(red=1,green=seq(from=1,to=0,length.out=50),blue=seq(from=1,to=0,length.out=50))
  #plot(1:50,1:50,col=colbarQuad1,pch=20,cex=2)
  colbarQuad2<-rgb(red=seq(from=1,to=0,length.out=50),green=0,blue=0)
  #plot(1:50,1:50,col=colbarQuad2,pch=20,cex=2)
  colbarQuad3<-rgb(red=0,green=0,blue=seq(from=0,to=1,length.out=50))
  #plot(1:50,1:50,col=colbarQuad3,pch=20,cex=2)
  colbarQuad4<-rgb(red=seq(from=0,to=1,length.out=50),green=seq(from=0,to=1,length.out=50),blue=1)
  #plot(1:50,1:50,col=colbarQuad4,pch=20,cex=2)
  colbar<-c(colbarQuad3,colbarQuad4,colbarQuad1,colbarQuad2)
  #plot(1:200,1:200,col=colbar,pch=20,cex=2)

  #***panels for fP1 
  htind<-5
  
  #short timescales
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0))
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  y<-Mod(fP1[inds])
  if (is.na(fP1_short_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP1_short_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(fP1[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)] 
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_short,xaxs="i",xaxt="n")
  mtext(latex2exp::TeX("$f_{P^{(1)}}$"),2,1.7,cex=2) 
  graphics::axis(1,at=log2(c(2,4,8)),labels=FALSE)
  text(xlimits_short[1],ylimits[2],ExplainInteractingMoran_PanLabs[1],adj=c(-0.1,1))
  
  #make the colorbar, which applies to all panels
  bds<-seq(from=.6*xlimits_short[1]+.4*xlimits_short[2],to=.1*xlimits_short[1]+.9*xlimits_short[2],
           length.out=length(breaks))
  rect(xleft=bds[1:(length(bds)-1)],
       ybottom=.1*ylimits[1]+.9*ylimits[2],
       xright=bds[2:length(bds)],
       ytop=ylimits[2],
       col=colbar,
       border=NA)
  ytags<-c(.1*ylimits[1]+.9*ylimits[2],.15*ylimits[1]+.85*ylimits[2],.175*ylimits[1]+.825*ylimits[2])
  lines(rep(.6*xlimits_short[1]+.4*xlimits_short[2],2),ytags[1:2])
  text(.6*xlimits_short[1]+.4*xlimits_short[2],ytags[3],latex2exp::TeX("$-\\pi$"),adj=c(.5,1))
  lines(rep((.6-.25)*xlimits_short[1]+(.4+.25)*xlimits_short[2],2),ytags[1:2])
  text((.6-.25)*xlimits_short[1]+(.4+.25)*xlimits_short[2],ytags[3],"0",adj=c(.5,1))
  lines(rep(.1*xlimits_short[1]+.9*xlimits_short[2],2),ytags[1:2])
  text(.1*xlimits_short[1]+.9*xlimits_short[2],ytags[3],latex2exp::TeX("$\\pi$"),adj=c(.5,1))
  
  #mid timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  y<-Mod(fP1[inds])
  if (is.na(fP1_mid_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP1_mid_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(fP1[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_mid,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(8,16)),labels=FALSE)
  text(xlimits_mid[1],ylimits[2],ExplainInteractingMoran_PanLabs[6],adj=c(-0.1,1))
  
  #long timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  y<-Mod(fP1[inds])
  if (is.na(fP1_long_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP1_long_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(fP1[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_long,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=FALSE)
  text(xlimits_long[1],ylimits[2],ExplainInteractingMoran_PanLabs[11],adj=c(-0.1,1))
  
  #***panels for fP2 
  htind<-4
  
  #short timescales
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  y<-Mod(fP2[inds])
  if (is.na(fP2_short_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP2_short_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(fP2[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)] 
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_short,xaxs="i",xaxt="n")
  mtext(latex2exp::TeX("$f_{P^{(2)}}$"),2,1.7,cex=2) 
  graphics::axis(1,at=log2(c(2,4,8)),labels=FALSE)
  text(xlimits_short[1],ylimits[2],ExplainInteractingMoran_PanLabs[2],adj=c(-0.1,1))
  
  #mid timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  y<-Mod(fP2[inds])
  if (is.na(fP2_mid_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP2_mid_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(fP2[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_mid,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(8,16)),labels=FALSE)
  text(xlimits_mid[1],ylimits[2],ExplainInteractingMoran_PanLabs[7],adj=c(-0.1,1))
  
  #long timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  y<-Mod(fP2[inds])
  if (is.na(fP2_long_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP2_long_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(fP2[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_long,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=FALSE)
  text(xlimits_long[1],ylimits[2],ExplainInteractingMoran_PanLabs[12],adj=c(-0.1,1))
  
  #***panels for fP1*Conj(fP2)
  htind<-3
  fP1ConjfP2<-fP1*Conj(fP2)
  
  #short timescales
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  y<-Mod(fP1ConjfP2[inds])
  if (is.na(fP1ConjfP2_short_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP1ConjfP2_short_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(fP1ConjfP2[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)] 
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_short,xaxs="i",xaxt="n")
  mtext(latex2exp::TeX("$f_{P^{(1)}} \\bar{f_{P^{(2)}}}$"),2,1.7,cex=2) 
  graphics::axis(1,at=log2(c(2,4,8)),labels=FALSE)
  text(xlimits_short[1],ylimits[2],ExplainInteractingMoran_PanLabs[3],adj=c(-0.1,1))
  
  #mid timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  y<-Mod(fP1ConjfP2[inds])
  if (is.na(fP1ConjfP2_mid_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP1ConjfP2_mid_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(fP1ConjfP2[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_mid,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(8,16)),labels=FALSE)
  text(xlimits_mid[1],ylimits[2],ExplainInteractingMoran_PanLabs[8],adj=c(-0.1,1))
  
  #long timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  y<-Mod(fP1ConjfP2[inds])
  if (is.na(fP1ConjfP2_long_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP1ConjfP2_long_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(fP1ConjfP2[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_long,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=FALSE)
  text(xlimits_long[1],ylimits[2],ExplainInteractingMoran_PanLabs[13],adj=c(-0.1,1))
  
  #***panels for Se1e2_avg
  htind<-2

  #short timescales
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  y<-Mod(Se1e2_avg[inds])
  if (is.na(Se1e2_short_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Se1e2_short_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(Se1e2_avg[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)] 
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_short,xaxs="i",xaxt="n")
  mtext(latex2exp::TeX("$\\rho_{\\epsilon^{(1)}\\epsilon^{(2)}}$"),2,1.7,cex=2) 
  graphics::axis(1,at=log2(c(2,4,8)),labels=FALSE)
  text(xlimits_short[1],ylimits[2],ExplainInteractingMoran_PanLabs[4],adj=c(-0.1,1))
  
  #mid timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  y<-Mod(Se1e2_avg[inds])
  if (is.na(Se1e2_mid_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Se1e2_mid_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(Se1e2_avg[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_mid,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(8,16)),labels=FALSE)
  text(xlimits_mid[1],ylimits[2],ExplainInteractingMoran_PanLabs[9],adj=c(-0.1,1))
  
  #long timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  y<-Mod(Se1e2_avg[inds])
  if (is.na(Se1e2_long_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Se1e2_long_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(Se1e2_avg[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_long,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=FALSE)
  text(xlimits_long[1],ylimits[2],ExplainInteractingMoran_PanLabs[14],adj=c(-0.1,1))
  
  #***panels for fP1*Conj(fP2)*Se1e2_avg
  htind<-1
  Prod<-fP1*Conj(fP2)*Se1e2_avg
  
  #short timescales
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  y<-Mod(Prod[inds])
  if (is.na(Prod_short_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Prod_short_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(Prod[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)] 
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_short,xaxs="i",xaxt="n")
  mtext(latex2exp::TeX("$f_{P^{(1)}} \\bar{f_{P^{(2)}}} \\rho_{\\epsilon^{(1)}\\epsilon^{(2)}}$"),2,1.7,cex=2) 
  graphics::axis(1,at=log2(c(2,4,8)),labels=c(0.5,1,2))
  text(xlimits_short[1],ylimits[2],ExplainInteractingMoran_PanLabs[5],adj=c(-0.1,1))
  
  #mid timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  y<-Mod(Prod[inds])
  if (is.na(Prod_mid_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Prod_mid_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(Prod[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_mid,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(8,16)),labels=c(2,4))
  mtext("Timescale, years",1,2.3,cex=2)
  text(xlimits_mid[1],ylimits[2],ExplainInteractingMoran_PanLabs[10],adj=c(-0.1,1))
  
  #long timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  y<-Mod(Prod[inds])
  if (is.na(Prod_long_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Prod_long_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  h<-Arg(Prod[inds])
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales[inds],y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xlim=xlimits_long,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=c(4,8,16,32))
  text(xlimits_long[1],ylimits[2],ExplainInteractingMoran_PanLabs[15],adj=c(-0.1,1))
  
  dev.off()
}

#A function for making the manuscript plots explaining direct Moran effects
#
explain_direct_Moran_ForMS<-function(freq,frg,fP1,fP2,fB,Se1e1_avg,Se2e2_avg,fnpre,plottype,PlotForMSArgs)
{
  #unpack plotting params
  fP1_short_ylims<-PlotForMSArgs$ExplainDirectMoran_fP1_short_ylims
  fP1_mid_ylims<-PlotForMSArgs$ExplainDirectMoran_fP1_mid_ylims
  fP1_long_ylims<-PlotForMSArgs$ExplainDirectMoran_fP1_long_ylims
  Se1e1_short_ylims<-PlotForMSArgs$ExplainDirectMoran_Se1e1_short_ylims
  Se1e1_mid_ylims<-PlotForMSArgs$ExplainDirectMoran_Se1e1_mid_ylims
  Se1e1_long_ylims<-PlotForMSArgs$ExplainDirectMoran_Se1e1_long_ylims
  fP2_short_ylims<-PlotForMSArgs$ExplainDirectMoran_fP2_short_ylims
  fP2_mid_ylims<-PlotForMSArgs$ExplainDirectMoran_fP2_mid_ylims
  fP2_long_ylims<-PlotForMSArgs$ExplainDirectMoran_fP2_long_ylims
  Se2e2_short_ylims<-PlotForMSArgs$ExplainDirectMoran_Se2e2_short_ylims
  Se2e2_mid_ylims<-PlotForMSArgs$ExplainDirectMoran_Se2e2_mid_ylims
  Se2e2_long_ylims<-PlotForMSArgs$ExplainDirectMoran_Se2e2_long_ylims
  fB_short_ylims<-PlotForMSArgs$ExplainDirectMoran_fB_short_ylims
  fB_mid_ylims<-PlotForMSArgs$ExplainDirectMoran_fB_mid_ylims
  fB_long_ylims<-PlotForMSArgs$ExplainDirectMoran_fB_long_ylims
  ExplainDirectMoran_PanLabs<-PlotForMSArgs$ExplainDirectMoran_PanLabs
  
  #constraint to certain frequencies, for plotting
  fP1<-fP1[frg[1]<freq & freq<=frg[2]]
  fP2<-fP2[frg[1]<freq & freq<=frg[2]]
  fB<-fB[frg[1]<freq & freq<=frg[2]]
  Se1e1_avg<-Se1e1_avg[frg[1]<freq & freq<=frg[2]]
  Se2e2_avg<-Se2e2_avg[frg[1]<freq & freq<=frg[2]]
  freq<-freq[frg[1]<freq & freq<=frg[2]]
  
  #we are going to plot against log2 timescale, on three panels for different timescale ranges, so
  #prep for that
  timescales<-1/freq
  l2timescales<-log2(timescales)
  xlimits<-range(l2timescales)
  xlimits_short<-c(xlimits[1],log2(8))
  xlimits_mid<-c(log2(8),log2(16))
  xlimits_long<-c(log2(16),xlimits[2])
  
  #panel layout quantities
  gap<-.15
  yaxnumwd<-.4
  yaxlabwd<-.4
  yaxwd<-yaxnumwd+yaxlabwd
  xaxht<-yaxwd
  panht<-1.25
  panwd_short<-2
  panwd_mid<-1
  panwd_long<-2
  totpanwd<-panwd_short+panwd_mid+panwd_long
  totht<-xaxht+5*panht+5*gap
  totwd<-yaxwd+2*yaxnumwd+totpanwd+gap
  
  #get the plotting device
  if (plottype=="jpg")
  {
    jpeg(paste0(resloc,fnpre,"_ExplainDirectMoran.jpg"),quality=95,width=totwd,height=totht,units="in",res=300)
  }
  if (plottype=="pdf")
  {
    pdf(paste0(resloc,fnpre,"_ExplainDirectMoran.pdf"),width=totwd,height=totht)
  }
  
  #***panels for fP1 
  htind<-5
  
  #short timescales
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0))
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  y<-(Mod(fP1[inds]))^2
  if (is.na(fP1_short_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP1_short_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_short,xaxs="i",xaxt="n")
  mtext(latex2exp::TeX("$|f_{P^{(1)}}|^2$"),2,1.7,cex=2) 
  graphics::axis(1,at=log2(c(2,4,8)),labels=FALSE)
  text(xlimits_short[1],ylimits[2],ExplainDirectMoran_PanLabs[1],adj=c(-0.1,1))
  
  #mid timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  y<-(Mod(fP1[inds]))^2
  if (is.na(fP1_mid_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP1_mid_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_mid,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(8,16)),labels=FALSE)
  text(xlimits_mid[1],ylimits[2],ExplainDirectMoran_PanLabs[6],adj=c(-0.1,1))
  
  #long timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  y<-(Mod(fP1[inds]))^2
  if (is.na(fP1_long_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP1_long_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_long,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=FALSE)
  text(xlimits_long[1],ylimits[2],ExplainDirectMoran_PanLabs[11],adj=c(-0.1,1))
  
  #***panels for Se1e1
  htind<-4
    
  #short timescales
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  y<-Se1e1_avg[inds]
  if (is.na(Se1e1_short_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Se1e1_short_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_short,xaxs="i",xaxt="n")
  mtext(latex2exp::TeX("$\\rho_{\\epsilon^{(1)}\\epsilon^{(1)}}$"),2,1.7,cex=2) 
  graphics::axis(1,at=log2(c(2,4,8)),labels=FALSE)
  text(xlimits_short[],ylimits[2],ExplainDirectMoran_PanLabs[2],adj=c(-0.1,1))
  
  #mid timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  y<-Se1e1_avg[inds]
  if (is.na(Se1e1_mid_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Se1e1_mid_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_mid,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(8,16)),labels=FALSE)
  text(xlimits_mid[1],ylimits[2],ExplainDirectMoran_PanLabs[7],adj=c(-0.1,1))
  
  #long timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  y<-Se1e1_avg[inds]
  if (is.na(Se1e1_long_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Se1e1_long_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_long,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=FALSE)
  text(xlimits_long[1],ylimits[2],ExplainDirectMoran_PanLabs[12],adj=c(-0.1,1))
  
  #***panels for fP2
  htind<-3
  
  #short timescales
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  y<-(Mod(fP2[inds]))^2
  if (is.na(fP2_short_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP2_short_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_short,xaxs="i",xaxt="n")
  mtext(latex2exp::TeX("$|f_{P^{(2)}}|^2$"),2,1.7,cex=2) 
  graphics::axis(1,at=log2(c(2,4,8)),labels=FALSE)
  text(xlimits_short[1],ylimits[2],ExplainDirectMoran_PanLabs[3],adj=c(-0.1,1))
  
  #mid timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  y<-(Mod(fP2[inds]))^2
  if (is.na(fP2_mid_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP2_mid_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_mid,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(8,16)),labels=FALSE)
  text(xlimits_mid[1],ylimits[2],ExplainDirectMoran_PanLabs[8],adj=c(-0.1,1))
  
  #long timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  y<-(Mod(fP2[inds]))^2
  if (is.na(fP2_long_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fP2_long_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_long,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=FALSE)
  text(xlimits_long[1],ylimits[2],ExplainDirectMoran_PanLabs[13],adj=c(-0.1,1))
  
  #***panels for Se2e2
  htind<-2
  
  #short timescales
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  y<-Se2e2_avg[inds]
  if (is.na(Se2e2_short_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Se2e2_short_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_short,xaxs="i",xaxt="n")
  mtext(latex2exp::TeX("$\\rho_{\\epsilon^{(2)}\\epsilon^{(2)}}$"),2,1.7,cex=2) 
  graphics::axis(1,at=log2(c(2,4,8)),labels=FALSE)
  text(xlimits_short[1],ylimits[2],ExplainDirectMoran_PanLabs[4],adj=c(-0.1,1))
  
  #mid timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  y<-Se2e2_avg[inds]
  if (is.na(Se2e2_mid_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Se2e2_mid_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_mid,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(8,16)),labels=FALSE)
  text(xlimits_mid[1],ylimits[2],ExplainDirectMoran_PanLabs[9],adj=c(-0.1,1))
  
  #long timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  y<-Se2e2_avg[inds]
  if (is.na(Se2e2_long_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-Se2e2_long_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_long,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=FALSE)
  text(xlimits_long[1],ylimits[2],ExplainDirectMoran_PanLabs[14],adj=c(-0.1,1))
  
  #***panels for fB
  htind<-1
  
  #short timescales
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  y<-1/((Mod(fB[inds]))^2)
  if (is.na(fB_short_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fB_short_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_short,xaxs="i",xaxt="n")
  mtext(latex2exp::TeX("$1/|f_B|^2$"),2,1.7,cex=2) 
  graphics::axis(1,at=log2(c(2,4,8)),labels=c(0.5,1,2))
  text(xlimits_short[1],ylimits[2],ExplainDirectMoran_PanLabs[5],adj=c(-0.1,1))
  
  #mid timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  y<-1/((Mod(fB[inds]))^2)
  if (is.na(fB_mid_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fB_mid_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_mid,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(8,16)),labels=c(2,4))
  mtext("Timescale, years",1,2.3,cex=2)
  text(xlimits_mid[1],ylimits[2],ExplainDirectMoran_PanLabs[10],adj=c(-0.1,1))
  
  #long timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht+(htind-1)*(panht+gap))/totht,
            (xaxht+(htind-1)*(panht+gap)+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  y<-1/((Mod(fB[inds]))^2)
  if (is.na(fB_long_ylims[1]))
  {
    ylimits<-range(y)
  }else
  {
    ylimits<-fB_long_ylims    
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],y,type="l",col="black",ylim=ylimits,xlim=xlimits_long,xaxs="i",xaxt="n")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=c(4,8,16,32))
  text(xlimits_long[1],ylimits[2],ExplainDirectMoran_PanLabs[15],adj=c(-0.1,1))
  
  dev.off()
}

#Utility function called in the above function, do_analysis
#
get_stats<-function(freq,T1_avg,T2_avg,T4_avg,totsync)
{
  #annual timescale stuff
  h<-c()
  inds<-which(freq>.2 & freq<=.3)
  h[1]<-mean(T1_avg[inds])  
  h[2]<-mean(T2_avg[inds])
  h[3]<-mean(T4_avg[inds])
  h[4]<-mean(totsync[inds])
  names(h)<-c("Ann sync NO3","Ann sync waves","Ann sync NO3/waves","Ann sync tot")
  res<-h
  
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
 
  return(res) 
}

#A utility function, called in the above function, do_analysis
#
make_plots_noise_sync<-function(freq,frg,S_avg,fnpre,plottype)
{
  #constraint to certain frequencies, for plotting
  S_avg<-S_avg[frg[1]<freq & freq<=frg[2]]
  freq<-freq[frg[1]<freq & freq<=frg[2]]
  
  #we are going to plot against log2 timescale, on three panels for different timescale ranges, so
  #prep for that
  timescales<-1/freq
  l2timescales<-log2(timescales)
  xlimits<-range(l2timescales)
  xlimits_short<-c(xlimits[1],log2(8))
  xlimits_mid<-c(log2(8),log2(16))
  xlimits_long<-c(log2(16),xlimits[2])
  
  #panel layout quantities
  gap<-.15
  yaxnumwd<-.4
  yaxlabwd<-.2
  yaxwd<-yaxnumwd+yaxlabwd
  xaxht<-yaxwd
  panht<-3
  panwd_short<-2
  panwd_mid<-1
  panwd_long<-2
  totpanwd<-panwd_short+panwd_mid+panwd_long
  totht<-xaxht+panht+gap
  totwd<-yaxwd+2*yaxnumwd+totpanwd+gap
  
  #get the plotting device for NO3
  if (plottype=="jpg")
  {
    jpeg(paste0(resloc,fnpre,"SyncPlot.jpg"),quality=95,width=totwd,height=totht,units="in",res=300)
  }
  if (plottype=="pdf")
  {
    pdf(paste0(resloc,fnpre,"SyncPlot.pdf"),width=totwd,height=totht)
  }
  
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0))
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  ylimits<-range(S_avg[inds])
  plot(l2timescales[inds],S_avg[inds],type="l",xaxt="n",col="black",
       ylim=ylimits,xlim=xlimits_short,xaxs="i")
  mtext("Synchrony",2,1.7)
  graphics::axis(1,at=log2(c(2,4,8)),labels=c(0.5,1,2))
  lines(l2timescales[inds],rep(0,length(timescales[inds])))
  lines(rep(log2(4),2),ylimits) #plot a vertical line at the annual frequency
  
  #the panel for mid (2-4 yr) timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  ylimits<-range(S_avg[inds])
  plot(l2timescales[inds],S_avg[inds],type="l",xaxt="n",col="black",
       ylim=ylimits,xlim=xlimits_mid,xaxs="i")
  mtext("Timescale, years",1,1.7)
  graphics::axis(1,at=log2(c(8,16)),labels=c(2,4))
  lines(l2timescales[inds],rep(0,length(timescales[inds])))
  
  #the panel for long (>4 yr) timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  ylimits<-range(S_avg[inds])
  plot(l2timescales[inds],S_avg[inds],type="l",xaxt="n",col="black",
       ylim=ylimits,xlim=xlimits_long,xaxs="i")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=c(4,8,16,32))
  lines(l2timescales[inds],rep(0,length(timescales[inds])))
  
  dev.off()
  
  return(NULL)
}

#A utility function, called in the above function, do_analysis. This is for making the plot that shows
#all the different contributors to synchrony.
#
make_plots<-function(freq,frg,totsync,T1_avg,T2_avg,T3_avg,T4_avg,T5_avg,T6_avg,resloc,fnpre,plottype)
{
  #constraint to certain frequencies, for plotting
  totsync<-totsync[frg[1]<freq & freq<=frg[2]]
  T1_avg<-T1_avg[frg[1]<freq & freq<=frg[2]]
  T2_avg<-T2_avg[frg[1]<freq & freq<=frg[2]]
  T3_avg<-T3_avg[frg[1]<freq & freq<=frg[2]]
  T4_avg<-T4_avg[frg[1]<freq & freq<=frg[2]]
  T5_avg<-T5_avg[frg[1]<freq & freq<=frg[2]]
  T6_avg<-T6_avg[frg[1]<freq & freq<=frg[2]]
  freq<-freq[frg[1]<freq & freq<=frg[2]]
  
  #create some other quantities to plot
  sumexpl<-T1_avg+T2_avg+T4_avg
  sumunexpl<-T3_avg+T5_avg+T6_avg
  sumterms<-sumexpl+sumunexpl
  
  #we are going to plot against log2 timescale, on three panels for different timescale ranges, so
  #prep for that
  timescales<-1/freq
  l2timescales<-log2(timescales)
  xlimits<-range(l2timescales)
  xlimits_short<-c(xlimits[1],log2(8))
  xlimits_mid<-c(log2(8),log2(16))
  xlimits_long<-c(log2(16),xlimits[2])
  
  #panel layout quantities
  gap<-.15
  yaxnumwd<-.4
  yaxlabwd<-.2
  yaxwd<-yaxnumwd+yaxlabwd
  xaxht<-yaxwd
  panht<-3
  panwd_short<-2
  panwd_mid<-1
  panwd_long<-2
  totpanwd<-panwd_short+panwd_mid+panwd_long
  totht<-xaxht+panht+gap
  totwd<-yaxwd+2*yaxnumwd+totpanwd+gap
  
  #get the plotting device
  if (plottype=="jpg")
  {
    jpeg(paste0(resloc,fnpre,"_MainPlot.jpg"),quality=95,width=totwd,height=totht,units="in",res=300)
  }
  if (plottype=="pdf")
  {
    pdf(paste0(resloc,fnpre,"_MainPlot.pdf"),width=totwd,height=totht)
  }
  
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0))
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  ylimits<-range(totsync[inds],T1_avg[inds],T2_avg[inds],T4_avg[inds],sumexpl[inds],sumterms[inds])
  plot(l2timescales[inds],totsync[inds],type="l",xaxt="n",col="black",
       ylim=ylimits,xlim=xlimits_short,xaxs="i")
  mtext("Synchrony contribution",2,1.7)
  graphics::axis(1,at=log2(c(2,4,8)),labels=c(0.5,1,2))
  lines(l2timescales[inds],T1_avg[inds],type="l",col="green") #direct effects of NO3
  lines(l2timescales[inds],T2_avg[inds],type="l",col="blue") #direct effects of waves
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="solid",col="green") 
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="dashed",col="blue") #interactions between waves and NO3
  lines(l2timescales[inds],sumexpl[inds],type="l",col="red",lty="solid") #explained stuff - sum of all the terms relating only to NO3 and waves (direct NO3 and wave effects, and interactions between those two)
  lines(l2timescales[inds],sumterms[inds],type="l",lty="dashed") #total of all terms, should approx equal the plot of totsync
  #lines(l2timescales[inds],sumunexpl[inds],type="l",lty="solid",col="red") #unexplained
  lines(l2timescales[inds],rep(0,length(timescales[inds])))
  lines(rep(log2(4),2),ylimits) #plot a vertical line at the annual frequency
  
  #the panel for mid (2-4 yr) timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  ylimits<-range(totsync[inds],T1_avg[inds],T2_avg[inds],T4_avg[inds],sumexpl[inds],sumterms[inds])
  plot(l2timescales[inds],totsync[inds],type="l",xaxt="n",col="black",
       ylim=ylimits,xlim=xlimits_mid,xaxs="i")
  mtext("Timescale, years",1,1.7)
  graphics::axis(1,at=log2(c(8,16)),labels=c(2,4))
  lines(l2timescales[inds],T1_avg[inds],type="l",col="green") #direct effects of NO3
  lines(l2timescales[inds],T2_avg[inds],type="l",col="blue") #direct effects of waves
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="solid",col="green") 
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="dashed",col="blue") #interactions between waves and NO3
  lines(l2timescales[inds],sumexpl[inds],type="l",col="red",lty="solid") #explained stuff - sum of all the terms relating only to NO3 and waves (direct NO3 and wave effects, and interactions between those two)
  lines(l2timescales[inds],sumterms[inds],type="l",lty="dashed") #total of all terms, should approx equal the plot of totsync
  #lines(l2timescales[inds],sumunexpl[inds],type="l",lty="solid",col="red") #unexplained
  lines(l2timescales[inds],rep(0,length(timescales[inds])))
  
  #the panel for long (>4 yr) timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  ylimits<-range(totsync[inds],T1_avg[inds],T2_avg[inds],T4_avg[inds],sumexpl[inds],sumterms[inds])
  plot(l2timescales[inds],totsync[inds],type="l",xaxt="n",col="black",
       ylim=ylimits,xlim=xlimits_long,xaxs="i")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=c(4,8,16,32))
  lines(l2timescales[inds],T1_avg[inds],type="l",col="green") #direct effects of NO3
  lines(l2timescales[inds],T2_avg[inds],type="l",col="blue") #direct effects of waves
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="solid",col="green") 
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="dashed",col="blue") #interactions between waves and NO3
  lines(l2timescales[inds],sumexpl[inds],type="l",col="red",lty="solid") #explained stuff - sum of all the terms relating only to NO3 and waves (direct NO3 and wave effects, and interactions between those two)
  lines(l2timescales[inds],sumterms[inds],type="l",lty="dashed") #total of all terms, should approx equal the plot of totsync
  #lines(l2timescales[inds],sumunexpl[inds],type="l",lty="solid",col="red") #unexplained
  lines(l2timescales[inds],rep(0,length(timescales[inds])))
  
  dev.off()
  
  return(NULL)
}

#A utility function, called in the above function, do_analysis. This is for making the plot that shows
#all the different contributors to synchrony, but presentation quality, for the manuscript.
#
make_plots_ForMS<-function(freq,frg,totsync,T1_avg,T2_avg,T3_avg,T4_avg,T5_avg,T6_avg,resloc,fnpre,plottype,PlotForMSArgs)
{
  #extract the presentational parameters you need 
  MainPlot_PanLabs<-PlotForMSArgs$MainPlot_PanLabs
  MainPlot_short_ylim<-PlotForMSArgs$MainPlot_short_ylim
  MainPlot_mid_ylim<-PlotForMSArgs$MainPlot_mid_ylim
  MainPlot_long_ylim<-PlotForMSArgs$MainPlot_long_ylim
    
  #constraint to certain frequencies, for plotting
  totsync<-totsync[frg[1]<freq & freq<=frg[2]]
  T1_avg<-T1_avg[frg[1]<freq & freq<=frg[2]]
  T2_avg<-T2_avg[frg[1]<freq & freq<=frg[2]]
  T3_avg<-T3_avg[frg[1]<freq & freq<=frg[2]]
  T4_avg<-T4_avg[frg[1]<freq & freq<=frg[2]]
  T5_avg<-T5_avg[frg[1]<freq & freq<=frg[2]]
  T6_avg<-T6_avg[frg[1]<freq & freq<=frg[2]]
  freq<-freq[frg[1]<freq & freq<=frg[2]]
  
  #create some other quantities to plot
  sumexpl<-T1_avg+T2_avg+T4_avg
  sumunexpl<-T3_avg+T5_avg+T6_avg
  sumterms<-sumexpl+sumunexpl
  
  #we are going to plot against log2 timescale, on three panels for different timescale ranges, so
  #prep for that
  timescales<-1/freq
  l2timescales<-log2(timescales)
  xlimits<-range(l2timescales)
  xlimits_short<-c(xlimits[1],log2(8))
  xlimits_mid<-c(log2(8),log2(16))
  xlimits_long<-c(log2(16),xlimits[2])
  
  #panel layout quantities
  gap<-.15
  yaxnumwd<-.4
  yaxlabwd<-.2
  yaxwd<-yaxnumwd+yaxlabwd
  xaxht<-yaxwd
  panht<-2
  panwd_short<-2
  panwd_mid<-1
  panwd_long<-2
  totpanwd<-panwd_short+panwd_mid+panwd_long
  totht<-xaxht+panht+gap
  totwd<-yaxwd+2*yaxnumwd+totpanwd+gap
  
  #get the plotting device
  if (plottype=="jpg")
  {
    jpeg(paste0(resloc,fnpre,"_MainPlot_ForMS.jpg"),quality=95,width=totwd,height=totht,units="in",res=300)
  }
  if (plottype=="pdf")
  {
    pdf(paste0(resloc,fnpre,"_MainPlot_ForMS.pdf"),width=totwd,height=totht)
  }
  
  #the panel for short (annual) timescales
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0))
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  if (is.na(MainPlot_short_ylim[1]))
  {
    ylimits<-range(totsync[inds],T1_avg[inds],T2_avg[inds],T4_avg[inds],sumexpl[inds],sumterms[inds])
  }else
  {
    ylimits<-MainPlot_short_ylim
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],totsync[inds],type="l",xaxt="n",col="black",
       ylim=ylimits,xlim=xlimits_short,xaxs="i")
  text(xlimits_short[1],ylimits[2],MainPlot_PanLabs[1],adj=c(-.1,1),cex=1.5)
  mtext("Synchrony contribution",2,1.7)
  graphics::axis(1,at=log2(c(2,4,8)),labels=c(0.5,1,2))
  lines(l2timescales[inds],T1_avg[inds],type="l",col="green") #direct effects of NO3
  lines(l2timescales[inds],T2_avg[inds],type="l",col="blue") #direct effects of waves
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="solid",col="green") 
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="dashed",col="blue") #interactions between waves and NO3
  lines(l2timescales[inds],sumexpl[inds],type="l",col="red",lty="solid") #explained stuff - sum of all the terms relating only to NO3 and waves (direct NO3 and wave effects, and interactions between those two)
  lines(l2timescales[inds],sumterms[inds],type="l",lty="dashed") #total of all terms, should approx equal the plot of totsync
  #lines(l2timescales[inds],sumunexpl[inds],type="l",lty="solid",col="red") #unexplained
  lines(l2timescales[inds],rep(0,length(timescales[inds])))
  lines(rep(log2(4),2),ylimits) #plot a vertical line at the annual frequency
  
  #the panel for mid (2-4 yr) timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  if (is.na(MainPlot_mid_ylim[1]))
  {
    ylimits<-range(totsync[inds],T1_avg[inds],T2_avg[inds],T4_avg[inds],sumexpl[inds],sumterms[inds])
  }else
  {
    ylimits<-MainPlot_mid_ylim
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],totsync[inds],type="l",xaxt="n",col="black",
       ylim=ylimits,xlim=xlimits_mid,xaxs="i")
  text(xlimits_mid[1],ylimits[2],MainPlot_PanLabs[2],adj=c(-.1,1),cex=1.5)
  mtext("Timescale, years",1,1.7)
  graphics::axis(1,at=log2(c(8,16)),labels=c(2,4))
  lines(l2timescales[inds],T1_avg[inds],type="l",col="green") #direct effects of NO3
  lines(l2timescales[inds],T2_avg[inds],type="l",col="blue") #direct effects of waves
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="solid",col="green") 
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="dashed",col="blue") #interactions between waves and NO3
  lines(l2timescales[inds],sumexpl[inds],type="l",col="red",lty="solid") #explained stuff - sum of all the terms relating only to NO3 and waves (direct NO3 and wave effects, and interactions between those two)
  lines(l2timescales[inds],sumterms[inds],type="l",lty="dashed") #total of all terms, should approx equal the plot of totsync
  #lines(l2timescales[inds],sumunexpl[inds],type="l",lty="solid",col="red") #unexplained
  lines(l2timescales[inds],rep(0,length(timescales[inds])))
  
  #the panel for long (>4 yr) timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  if (is.na(MainPlot_long_ylim[1]))
  {
    ylimits<-range(totsync[inds],T1_avg[inds],T2_avg[inds],T4_avg[inds],sumexpl[inds],sumterms[inds])
  }else
  {
    ylimits<-MainPlot_long_ylim
  }
  ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  plot(l2timescales[inds],totsync[inds],type="l",xaxt="n",col="black",
       ylim=ylimits,xlim=xlimits_long,xaxs="i")
  text(xlimits_long[1],ylimits[2],MainPlot_PanLabs[3],adj=c(-.1,1),cex=1.5)
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=c(4,8,16,32))
  lines(l2timescales[inds],T1_avg[inds],type="l",col="green") #direct effects of NO3
  lines(l2timescales[inds],T2_avg[inds],type="l",col="blue") #direct effects of waves
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="solid",col="green") 
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="dashed",col="blue") #interactions between waves and NO3
  lines(l2timescales[inds],sumexpl[inds],type="l",col="red",lty="solid") #explained stuff - sum of all the terms relating only to NO3 and waves (direct NO3 and wave effects, and interactions between those two)
  lines(l2timescales[inds],sumterms[inds],type="l",lty="dashed") #total of all terms, should approx equal the plot of totsync
  #lines(l2timescales[inds],sumunexpl[inds],type="l",lty="solid",col="red") #unexplained
  lines(l2timescales[inds],rep(0,length(timescales[inds])))
  
  dev.off()
  
  return(NULL)
}

#A utility function, called in the above function, do_analysis. Plots a complex vector's
#phase and magnitude separately.
#
make_plot_component<-function(freq,frg,comp,plotname,plottype)
{
  #constraint to certain frequencies, for plotting
  comp<-comp[frg[1]<freq & freq<=frg[2]]
  freq<-freq[frg[1]<freq & freq<=frg[2]]
  
  #we are going to plot against log2 timescale, on three panels for different timescale ranges, so
  #prep for that
  timescales<-1/freq
  l2timescales<-log2(timescales)
  xlimits<-range(l2timescales)
  xlimits_short<-c(xlimits[1],log2(8))
  xlimits_mid<-c(log2(8),log2(16))
  xlimits_long<-c(log2(16),xlimits[2])
  
  l2timescales<-rev(l2timescales)
  comp<-rev(comp)
  
  #panel layout quantities
  gap<-.25
  yaxnumwd<-.4
  yaxlabwd<-.2
  yaxwd<-yaxnumwd+yaxlabwd
  xaxht<-yaxwd
  panht<-1.5
  panwd_short<-2
  panwd_mid<-1
  panwd_long<-2
  totpanwd<-panwd_short+panwd_mid+panwd_long
  totht<-xaxht+2*panht+2*gap
  totwd<-yaxwd+2*yaxnumwd+totpanwd+gap
  
  #get the plotting device
  if (plottype=="jpg")
  {
    jpeg(paste0(resloc,plotname,".jpg"),quality=95,width=totwd,height=totht,units="in",res=300)
  }
  if (plottype=="pdf")
  {
    pdf(paste0(resloc,plotname,".pdf"),width=totwd,height=totht)
  }
  
  #the panel for short (annual) timescales, phase
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht+panht+gap)/totht,
            (xaxht+2*panht+gap)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0))
  inds<-which(l2timescales>=xlimits_short[1] & l2timescales<=xlimits_short[2])
  plot(0,0,xaxt="n",type="n",xlim=xlimits_short,ylim=c(-pi,pi),xaxs="i")
  plotphasefunc(l2timescales[inds],comp[inds],pch=20,cex=0.5,lty="solid",col="black")
  lines(xlimits_short,rep(pi/2,2),type="l",lty="dashed")
  lines(xlimits_short,rep(-pi/2,2),type="l",lty="dashed")
  mtext("Phase",2,1.7)
  graphics::axis(1,at=log2(c(2,4,8)),labels=FALSE)
  
  #the panel for short (annual) timescales, magnitude
  par(fig=c((yaxwd)/totwd,
            (yaxwd+panwd_short)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  plot(l2timescales[inds],Mod(comp[inds]),type="p",xaxt="n",
       xlim=xlimits_short,xaxs="i",pch=20,cex=0.5,col="black")
  lines(l2timescales[inds],Mod(comp[inds]),type="l",lty="solid",col="black")
  mtext("Modulus",2,1.7)
  graphics::axis(1,at=log2(c(2,4,8)),labels=c(0.5,1,2))
  
  #the panel for mid (2-4 yr) timescales, phase
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht+panht+gap)/totht,
            (xaxht+2*panht+gap)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_mid[1] & l2timescales<=xlimits_mid[2])
  plot(0,0,xaxt="n",type="n",xlim=xlimits_mid,ylim=c(-pi,pi),xaxs="i")
  plotphasefunc(l2timescales[inds],comp[inds],pch=20,cex=0.5,lty="solid",col="black")
  lines(xlimits_mid,rep(pi/2,2),type="l",lty="dashed")
  lines(xlimits_mid,rep(-pi/2,2),type="l",lty="dashed")
  graphics::axis(1,at=log2(c(8,16)),labels=FALSE)
  
  #the panel for mid (2-4 yr) timescales, magnitude
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  plot(l2timescales[inds],Mod(comp[inds]),type="b",xaxt="n",
       xlim=xlimits_mid,xaxs="i",pch=20,cex=0.5,col="black")
  lines(l2timescales[inds],Mod(comp[inds]),type="l",lty="solid",col="black")
  graphics::axis(1,at=log2(c(8,16)),labels=c(2,4))
  mtext("Timescale, years",1,1.7)
  
  #the panel for long (>4 yr) timescales, phase
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht+panht+gap)/totht,
            (xaxht+2*panht+gap)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  plot(0,0,xaxt="n",type="n",xlim=xlimits_long,ylim=c(-pi,pi),xaxs="i")
  plotphasefunc(l2timescales[inds],comp[inds],pch=20,cex=0.5,lty="solid",col="black")
  lines(xlimits_long,rep(pi/2,2),type="l",lty="dashed")
  lines(xlimits_long,rep(-pi/2,2),type="l",lty="dashed")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=FALSE)
  
  #the panel for long (>4 yr) timescales, magnitude
  par(fig=c((yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_mid+yaxnumwd+panwd_long)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  plot(l2timescales[inds],Mod(comp[inds]),type="b",xaxt="n",
       xlim=xlimits_long,xaxs="i",pch=20,cex=0.5,col="black")
  lines(l2timescales[inds],Mod(comp[inds]),type="l",lty="solid",col="black")
  graphics::axis(1,at=log2(c(16,32,64,128)),labels=c(4,8,16,32))
  
  dev.off()

  return(NULL)
}

# #***
# #some tests of the main function, based on comparison to analyses done elsewhere using seperate code
# #***
# 
# Args<-list(locstouse=25:220,lags=c(4,1,0),frg=c(0,.35),fnpre="Test1",plottype="pdf")
# do_analysis(Args)
# Args<-list(locstouse=25:220,lags=c(8,1,0),frg=c(0,.35),fnpre="Test2",plottype="pdf")
# do_analysis(Args)
#
# #***
# #Do a big run of the function for sets of sites up and down the coast
# #***
# 
# #assemble a list of argument lists for a call to mclapply, and do that call
# alldists<-ncf::gcdist(locs$Lon,locs$Lat)
# allargs<-list()
# allargsind<-1
# centerinds<-c()
# for (counter in 1:(dim(locs)[1]))
# {
#   locstouse<-which(alldists[counter,]<25)
#   if (length(locstouse)<20){ next }
#   centerinds<-c(centerinds,counter)
#   curargs1<-list(locstouse=locstouse,lags=c(4,1,0),frg=c(0,.5),fnpre=paste0("KelpLag_04_Center_",counter),plottype="jpg")
#   curargs2<-list(locstouse=locstouse,lags=c(8,1,0),frg=c(0,.5),fnpre=paste0("KelpLag_08_Center_",counter),plottype="jpg")
#   curargs3<-list(locstouse=locstouse,lags=c(12,1,0),frg=c(0,.5),fnpre=paste0("KelpLag_12_Center_",counter),plottype="jpg")
#   allargs[[allargsind]]<-curargs1
#   allargs[[allargsind+1]]<-curargs2
#   allargs[[allargsind+2]]<-curargs3
#   allargsind<-allargsind+3
# }
# allres<-parallel::mclapply(FUN=do_analysis,X=allargs,mc.cores=10)
# 
# saveRDS(allargs,file=paste0(resloc,"allargs.Rds"))
# saveRDS(allres,file=paste0(resloc,"allres.Rds"))
#
# #***
# #now organize summary results 
# #***
# 
# #put the results in a data frame for lag 4
# allres4<-matrix(NA,length(allres)/3,length(allres[[1]]))
# colnames(allres4)<-names(allres[[1]])
# for (counter in seq(from=1,by=3,to=length(allres)))
# {
#   allres4[(counter+2)/3,]<-allres[[counter]]
# }
# saveRDS(allres4,paste0(resloc,"allres4.Rds"))
# 
# #same for lag 8
# allres8<-matrix(NA,length(allres)/3,length(allres[[2]]))
# colnames(allres8)<-names(allres[[2]])
# for (counter in seq(from=1,by=3,to=length(allres)))
# {
#   allres8[(counter+2)/3,]<-allres[[counter+1]]
# }
# saveRDS(allres8,paste0(resloc,"allres8.Rds"))
# 
# #same for lag 12
# allres12<-matrix(NA,length(allres)/3,length(allres[[3]]))
# colnames(allres12)<-names(allres[[3]])
# for (counter in seq(from=1,by=3,to=length(allres)))
# {
#   allres12[(counter+2)/3,]<-allres[[counter+2]]
# }
# saveRDS(allres12,paste0(resloc,"allres12.Rds"))
# 
# #get variables we can use as x axes of the plots below
# mdlocind<-c() #will be median of indices of locations used in each run
# mnloclat<-c() #will be mean of latitudes of locations used
# mnloclon<-c() #will be mean of longitudes of locations used
# for (counter in seq(from=1,by=2,to=length(allargs)))
# {
#   mdlocind[(counter+1)/2]<-median(allargs[[counter]]$locstouse)
#   mnloclat[(counter+1)/2]<-mean(locs$Lat[allargs[[counter]]$locstouse])
#   mnloclon[(counter+1)/2]<-mean(locs$Lon[allargs[[counter]]$locstouse])
# }
# ctrloclat<-locs$Lat[centerinds]
# ctrloclon<-locs$Lon[centerinds]
# 
# #***
# #make plots which can (hopefully) be interpreted to get a sense for whether results are more or less
# #consistent across central CA, and more or less consistent across southern CA
# #***
# 
# put_loc_lines<-function()
# {
#   text(PtConcInd,ylimits[2],"Pt. Conc.",srt=90,adj=c(1,0))
#   lines(rep(PtConcInd,2),ylimits)
#   
#   text(MontBInd,ylimits[2],"Mont. Bay",srt=90,adj=c(1,0))
#   lines(rep(MontBInd,2),ylimits,lty="dashed")
#   
#   text(CarBInd,ylimits[2],"Car. Bay",srt=90,adj=c(1,0))
#   lines(rep(CarBInd,2),ylimits,lty="dashed")
#   
#   text(MorBInd,ylimits[2],"Mor. Bay",srt=90,adj=c(1,0))
#   lines(rep(MorBInd,2),ylimits,lty="dashed")
#   
#   text(SBInd,ylimits[2],"SB",srt=90,adj=c(1,0))
#   lines(rep(SBInd,2),ylimits,lty="dashed")
#     
#   lines(rep(OxInd,2),ylimits,lty="dashed")
#   text(OxInd,ylimits[2],"Oxnard",srt=90,adj=c(1,0))
#   
#   lines(rep(LAInd,2),ylimits,lty="dashed")
#   text(LAInd,ylimits[2],"LA",srt=90,adj=c(1,0))
# }
# 
# #***kelp lag 4 plots
# 
# #make plots about R^2 - how well does the ARMA capture the dynamics
# pdf(paste0(resloc,"LinearModel_Rsq_v_centerind_kelplag4.pdf"))
# ylimits<-range(allres4[,'model r sq'])
# plot(centerinds,allres4[,'model r sq'],type='b',xlab="Center index",ylab="Linear model R sq")
# put_loc_lines()
# dev.off()
# 
# #make plots about coefficients - kelp coefficients first
# pdf(paste0(resloc,"LinearModel_KelpCoefs_v_centerind_kelplag4.pdf"))
# ylimits<-range(allres4[,'coef_kelp_l1'],allres4[,'coef_kelp_l2'],allres4[,'coef_kelp_l2'],allres4[,'coef_kelp_l4'])
# plot(centerinds,allres4[,'coef_kelp_l1'],type='b',xlab="Center index",ylab="Kelp coefficient",pch=1,
#      ylim=ylimits)
# lines(centerinds,allres4[,'coef_kelp_l2'],type='b',pch=2)
# lines(centerinds,allres4[,'coef_kelp_l3'],type='b',pch=22)
# lines(centerinds,allres4[,'coef_kelp_l4'],type='b',pch=11)
# put_loc_lines()
# dev.off()
# 
# #now NO3 coefficients
# pdf(paste0(resloc,"LinearModel_NO3Coefs_v_centerind_kelplag4.pdf"))
# ylimits<-range(allres4[,'coef_NO3_l0'],allres4[,'coef_NO3_l1'])
# plot(centerinds,allres4[,'coef_NO3_l0'],type='b',xlab="Center index",ylab="NO3 coefficient",pch=1,
#      ylim=ylimits)
# lines(centerinds,allres4[,'coef_NO3_l1'],type='b',pch=2)
# put_loc_lines()
# dev.off()
# 
# #now waves coefficient
# pdf(paste0(resloc,"LinearModel_wavesCoefs_v_centerind_kelplag4.pdf"))
# ylimits<-range(allres4[,'coef_waves_l0'])
# plot(centerinds,allres4[,'coef_waves_l0'],type='b',xlab="Center index",ylab="Waves coefficient",pch=1,
#      ylim=ylimits)
# put_loc_lines()
# dev.off()
# 
# #now plot info about annual-timescale sycnhrony
# pdf(paste0(resloc,"Synchrony_AnnualTimescale_v_centerind_kelplag4_cov.pdf"))
# ylimits<-range(0,allres4[,'cov_Ann sync tot'],allres4[,'cov_Ann sync NO3'],allres4[,'cov_Ann sync waves'],allres4[,'cov_Ann sync NO3/waves'])
# plot(centerinds,allres4[,'cov_Ann sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, annual timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres4[,'cov_Ann sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres4[,'cov_Ann sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres4[,'cov_Ann sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres4[,'cov_Ann sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres4[,'cov_Ann sync NO3']+allres4[,'cov_Ann sync waves']+allres4[,'cov_Ann sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# pdf(paste0(resloc,"Synchrony_AnnualTimescale_v_centerind_kelplag4_cor.pdf"))
# ylimits<-range(0,allres4[,'cor_Ann sync tot'],allres4[,'cor_Ann sync NO3'],allres4[,'cor_Ann sync waves'],allres4[,'cor_Ann sync NO3/waves'])
# plot(centerinds,allres4[,'cor_Ann sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, annual timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres4[,'cor_Ann sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres4[,'cor_Ann sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres4[,'cor_Ann sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres4[,'cor_Ann sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres4[,'cor_Ann sync NO3']+allres4[,'cor_Ann sync waves']+allres4[,'cor_Ann sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# #***DAN: IDEA, maybe the reason why synchrony is so much less in SoCal is because waves are much
# #less synchronous, *and* their effects are no longer aligned with the effects of NO3. That means
# #you lose both the wave and NO3/wave interaction effects on synch. 
# 
# #now do 2-4 year timescales
# pdf(paste0(resloc,"Synchrony_2to4yrTimescales_v_centerind_kelplag4_cov.pdf"))
# ylimits<-range(0,allres4[,'cov_2to4yr sync tot'],allres4[,'cov_2to4yr sync NO3'],allres4[,'cov_2to4yr sync waves'],allres4[,'cov_2to4yr sync NO3/waves'])
# plot(centerinds,allres4[,'cov_2to4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, 2-4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres4[,'cov_2to4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres4[,'cov_2to4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres4[,'cov_2to4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres4[,'cov_2to4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres4[,'cov_2to4yr sync NO3']+allres4[,'cov_2to4yr sync waves']+allres4[,'cov_2to4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# pdf(paste0(resloc,"Synchrony_2to4yrTimescales_v_centerind_kelplag4_cor.pdf"))
# ylimits<-range(0,allres4[,'cor_2to4yr sync tot'],allres4[,'cor_2to4yr sync NO3'],allres4[,'cor_2to4yr sync waves'],allres4[,'cor_2to4yr sync NO3/waves'])
# plot(centerinds,allres4[,'cor_2to4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, 2-4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres4[,'cor_2to4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres4[,'cor_2to4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres4[,'cor_2to4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres4[,'cor_2to4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres4[,'cor_2to4yr sync NO3']+allres4[,'cor_2to4yr sync waves']+allres4[,'cor_2to4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# #now do >4 yr timescales
# pdf(paste0(resloc,"Synchrony_gt4yrTimescales_v_centerind_kelplag4_cov.pdf"))
# ylimits<-range(0,allres4[,'cov_>4yr sync tot'],allres4[,'cov_>4yr sync NO3'],allres4[,'cov_>4yr sync waves'],allres4[,'cov_>4yr sync NO3/waves'])
# plot(centerinds,allres4[,'cov_>4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, >4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres4[,'cov_>4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres4[,'cov_>4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres4[,'cov_>4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres4[,'cov_>4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres4[,'cov_>4yr sync NO3']+allres4[,'cov_>4yr sync waves']+allres4[,'cov_>4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# pdf(paste0(resloc,"Synchrony_gt4yrTimescales_v_centerind_kelplag4_cor.pdf"))
# ylimits<-range(0,allres4[,'cor_>4yr sync tot'],allres4[,'cor_>4yr sync NO3'],allres4[,'cor_>4yr sync waves'],allres4[,'cor_>4yr sync NO3/waves'])
# plot(centerinds,allres4[,'cor_>4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, >4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres4[,'cor_>4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres4[,'cor_>4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres4[,'cor_>4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres4[,'cor_>4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres4[,'cor_>4yr sync NO3']+allres4[,'cor_>4yr sync waves']+allres4[,'cor_>4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# #plot numbers of sites used 
# pdf(paste0(resloc,"Numsitesused_v_centerind_kelplag4.pdf"))
# ylimits<-range(allres4[,'num locs'])
# plot(centerinds,allres4[,'num locs'],type='b',xlab="Center index",ylab="Number of sites used")
# lines(rep(PtConcInd,2),ylimits)
# lines(rep(SBInd,2),ylimits,lty="dashed")
# put_loc_lines()
# dev.off()
# 
# #***DAN: THOUGHTS: I think if I choose a region in Central CA and one in SoCal (centered slightly west 
# #of SB) and contrast them, that will be good. I think generally results are roughly consistent across 
# #Central CA and across SoCal from Pt Conception to Oxnard (and after that the sites are too spaced out 
# #to do this analsis). The idea would be to present one analysis from Central Cal and one from So Cal 
# #and then say those are representative of if we had used different areas. Maybe do two from Central CA, 
# #since it's bigger than the usable So Cal area. Then you would not have to say your selected area is
# #representative, since you'd be using the bulk of the sites. So split up Carmel Bay to Morro Bay into 
# #two sets of sites, and use Pt Conc to Oxnard for the So Cal sites.
# 
# #***kelp lag 8 plots
# 
# #make plots about R^2 - how well does the ARMA capture the dynamics
# pdf(paste0(resloc,"LinearModel_Rsq_v_centerind_kelplag8.pdf"))
# ylimits<-range(allres8[,'model r sq'])
# plot(centerinds,allres8[,'model r sq'],type='b',xlab="Center index",ylab="Linear model R sq")
# put_loc_lines()
# dev.off()
# 
# #make plots about coefficients - kelp coefficients first - don't bother, too many lags
# 
# #now NO3 coefficients
# pdf(paste0(resloc,"LinearModel_NO3Coefs_v_centerind_kelplag8.pdf"))
# ylimits<-range(allres8[,'coef_NO3_l0'],allres8[,'coef_NO3_l1'])
# plot(centerinds,allres8[,'coef_NO3_l0'],type='b',xlab="Center index",ylab="NO3 coefficient",pch=1,
#      ylim=ylimits)
# lines(centerinds,allres8[,'coef_NO3_l1'],type='b',pch=2)
# put_loc_lines()
# dev.off()
# 
# #now waves coefficient
# pdf(paste0(resloc,"LinearModel_wavesCoefs_v_centerind_kelplag8.pdf"))
# ylimits<-range(allres8[,'coef_waves_l0'])
# plot(centerinds,allres8[,'coef_waves_l0'],type='b',xlab="Center index",ylab="Waves coefficient",pch=1,
#      ylim=ylimits)
# put_loc_lines()
# dev.off()
# 
# #now plot info about annual-timescale sycnhrony
# pdf(paste0(resloc,"Synchrony_AnnualTimescale_v_centerind_kelplag8_cov.pdf"))
# ylimits<-range(0,allres8[,'cov_Ann sync tot'],allres8[,'cov_Ann sync NO3'],allres8[,'cov_Ann sync waves'],allres8[,'cov_Ann sync NO3/waves'])
# plot(centerinds,allres8[,'cov_Ann sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, annual timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres8[,'cov_Ann sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres8[,'cov_Ann sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres8[,'cov_Ann sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres8[,'cov_Ann sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres8[,'cov_Ann sync NO3']+allres8[,'cov_Ann sync waves']+allres8[,'cov_Ann sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# pdf(paste0(resloc,"Synchrony_AnnualTimescale_v_centerind_kelplag8_cor.pdf"))
# ylimits<-range(0,allres8[,'cor_Ann sync tot'],allres8[,'cor_Ann sync NO3'],allres8[,'cor_Ann sync waves'],allres8[,'cor_Ann sync NO3/waves'])
# plot(centerinds,allres8[,'cor_Ann sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, annual timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres8[,'cor_Ann sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres8[,'cor_Ann sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres8[,'cor_Ann sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres8[,'cor_Ann sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres8[,'cor_Ann sync NO3']+allres8[,'cor_Ann sync waves']+allres8[,'cor_Ann sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# #now do 2-4 year timescales
# pdf(paste0(resloc,"Synchrony_2to4yrTimescales_v_centerind_kelplag8_cov.pdf"))
# ylimits<-range(0,allres8[,'cov_2to4yr sync tot'],allres8[,'cov_2to4yr sync NO3'],allres8[,'cov_2to4yr sync waves'],allres8[,'cov_2to4yr sync NO3/waves'])
# plot(centerinds,allres8[,'cov_2to4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, 2-4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres8[,'cov_2to4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres8[,'cov_2to4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres8[,'cov_2to4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres8[,'cov_2to4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres8[,'cov_2to4yr sync NO3']+allres8[,'cov_2to4yr sync waves']+allres8[,'cov_2to4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# pdf(paste0(resloc,"Synchrony_2to4yrTimescales_v_centerind_kelplag8_cor.pdf"))
# ylimits<-range(0,allres8[,'cor_2to4yr sync tot'],allres8[,'cor_2to4yr sync NO3'],allres8[,'cor_2to4yr sync waves'],allres8[,'cor_2to4yr sync NO3/waves'])
# plot(centerinds,allres8[,'cor_2to4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, 2-4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres8[,'cor_2to4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres8[,'cor_2to4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres8[,'cor_2to4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres8[,'cor_2to4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres8[,'cor_2to4yr sync NO3']+allres8[,'cor_2to4yr sync waves']+allres8[,'cor_2to4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# #now do >4 yr timescales
# pdf(paste0(resloc,"Synchrony_gt4yrTimescales_v_centerind_kelplag8_cov.pdf"))
# ylimits<-range(0,allres8[,'cov_>4yr sync tot'],allres8[,'cov_>4yr sync NO3'],allres8[,'cov_>4yr sync waves'],allres8[,'cov_>4yr sync NO3/waves'])
# plot(centerinds,allres8[,'cov_>4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, >4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres8[,'cov_>4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres8[,'cov_>4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres8[,'cov_>4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres8[,'cov_>4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres8[,'cov_>4yr sync NO3']+allres8[,'cov_>4yr sync waves']+allres8[,'cov_>4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# pdf(paste0(resloc,"Synchrony_gt4yrTimescales_v_centerind_kelplag8_cor.pdf"))
# ylimits<-range(0,allres8[,'cor_>4yr sync tot'],allres8[,'cor_>4yr sync NO3'],allres8[,'cor_>4yr sync waves'],allres8[,'cor_>4yr sync NO3/waves'])
# plot(centerinds,allres8[,'cor_>4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, >4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres8[,'cor_>4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres8[,'cor_>4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres8[,'cor_>4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres8[,'cor_>4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres8[,'cor_>4yr sync NO3']+allres8[,'cor_>4yr sync waves']+allres8[,'cor_>4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# #plot numbers of sites used - should be exactly the same as the corresponding kelp lag 4 plot
# pdf(paste0(resloc,"Numsitesused_v_centerind_kelplag8.pdf"))
# ylimits<-range(allres8[,'num locs'])
# plot(centerinds,allres8[,'num locs'],type='b',xlab="Center index",ylab="Number of sites used")
# lines(rep(PtConcInd,2),ylimits)
# lines(rep(SBInd,2),ylimits,lty="dashed")
# put_loc_lines()
# dev.off()
# 
# #***kelp lag 12 plots
# 
# #make plots about R^2 - how well does the ARMA capture the dynamics
# pdf(paste0(resloc,"LinearModel_Rsq_v_centerind_kelplag12.pdf"))
# ylimits<-range(allres12[,'model r sq'])
# plot(centerinds,allres12[,'model r sq'],type='b',xlab="Center index",ylab="Linear model R sq")
# put_loc_lines()
# dev.off()
# 
# #make plots about coefficients - kelp coefficients first - don't bother, too many lags
# 
# #now NO3 coefficients
# pdf(paste0(resloc,"LinearModel_NO3Coefs_v_centerind_kelplag12.pdf"))
# ylimits<-range(allres12[,'coef_NO3_l0'],allres12[,'coef_NO3_l1'])
# plot(centerinds,allres12[,'coef_NO3_l0'],type='b',xlab="Center index",ylab="NO3 coefficient",pch=1,
#      ylim=ylimits)
# lines(centerinds,allres12[,'coef_NO3_l1'],type='b',pch=2)
# put_loc_lines()
# dev.off()
# 
# #now waves coefficient
# pdf(paste0(resloc,"LinearModel_wavesCoefs_v_centerind_kelplag12.pdf"))
# ylimits<-range(allres12[,'coef_waves_l0'])
# plot(centerinds,allres12[,'coef_waves_l0'],type='b',xlab="Center index",ylab="Waves coefficient",pch=1,
#      ylim=ylimits)
# put_loc_lines()
# dev.off()
# 
# #now plot info about annual-timescale sycnhrony
# pdf(paste0(resloc,"Synchrony_AnnualTimescale_v_centerind_kelplag12_cov.pdf"))
# ylimits<-range(0,allres12[,'cov_Ann sync tot'],allres12[,'cov_Ann sync NO3'],allres12[,'cov_Ann sync waves'],allres12[,'cov_Ann sync NO3/waves'])
# plot(centerinds,allres12[,'cov_Ann sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, annual timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres12[,'cov_Ann sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres12[,'cov_Ann sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres12[,'cov_Ann sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres12[,'cov_Ann sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres12[,'cov_Ann sync NO3']+allres12[,'cov_Ann sync waves']+allres12[,'cov_Ann sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# pdf(paste0(resloc,"Synchrony_AnnualTimescale_v_centerind_kelplag12_cor.pdf"))
# ylimits<-range(0,allres12[,'cor_Ann sync tot'],allres12[,'cor_Ann sync NO3'],allres12[,'cor_Ann sync waves'],allres12[,'cor_Ann sync NO3/waves'])
# plot(centerinds,allres12[,'cor_Ann sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, annual timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres12[,'cor_Ann sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres12[,'cor_Ann sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres12[,'cor_Ann sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres12[,'cor_Ann sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres12[,'cor_Ann sync NO3']+allres12[,'cor_Ann sync waves']+allres12[,'cor_Ann sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# #***DAN: IDEA: Maybe I need to add to my linear modelling analysis an analysis that assesses
# #the appropriate number of kelp lags to use? 
# 
# #now do 2-4 year timescales
# pdf(paste0(resloc,"Synchrony_2to4yrTimescales_v_centerind_kelplag12_cov.pdf"))
# ylimits<-range(0,allres12[,'cov_2to4yr sync tot'],allres12[,'cov_2to4yr sync NO3'],allres12[,'cov_2to4yr sync waves'],allres12[,'cov_2to4yr sync NO3/waves'])
# plot(centerinds,allres12[,'cov_2to4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, 2-4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres12[,'cov_2to4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres12[,'cov_2to4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres12[,'cov_2to4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres12[,'cov_2to4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres12[,'cov_2to4yr sync NO3']+allres12[,'cov_2to4yr sync waves']+allres12[,'cov_2to4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# pdf(paste0(resloc,"Synchrony_2to4yrTimescales_v_centerind_kelplag12_cor.pdf"))
# ylimits<-range(0,allres12[,'cor_2to4yr sync tot'],allres12[,'cor_2to4yr sync NO3'],allres12[,'cor_2to4yr sync waves'],allres12[,'cor_2to4yr sync NO3/waves'])
# plot(centerinds,allres12[,'cor_2to4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, 2-4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres12[,'cor_2to4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres12[,'cor_2to4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres12[,'cor_2to4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres12[,'cor_2to4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres12[,'cor_2to4yr sync NO3']+allres12[,'cor_2to4yr sync waves']+allres12[,'cor_2to4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# #now do >4 yr timescales
# pdf(paste0(resloc,"Synchrony_gt4yrTimescales_v_centerind_kelplag12_cov.pdf"))
# ylimits<-range(0,allres12[,'cov_>4yr sync tot'],allres12[,'cov_>4yr sync NO3'],allres12[,'cov_>4yr sync waves'],allres12[,'cov_>4yr sync NO3/waves'])
# plot(centerinds,allres12[,'cov_>4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, >4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres12[,'cov_>4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres12[,'cov_>4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres12[,'cov_>4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres12[,'cov_>4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres12[,'cov_>4yr sync NO3']+allres12[,'cov_>4yr sync waves']+allres12[,'cov_>4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# pdf(paste0(resloc,"Synchrony_gt4yrTimescales_v_centerind_kelplag12_cor.pdf"))
# ylimits<-range(0,allres12[,'cor_>4yr sync tot'],allres12[,'cor_>4yr sync NO3'],allres12[,'cor_>4yr sync waves'],allres12[,'cor_>4yr sync NO3/waves'])
# plot(centerinds,allres12[,'cor_>4yr sync tot'],type='b',col="black",xlab="Center index",ylab="Component of synchrony, >4 yr timescales",pch=20,
#      ylim=ylimits)
# lines(centerinds,allres12[,'cor_>4yr sync NO3'],type='b',pch=20,col="green")
# lines(centerinds,allres12[,'cor_>4yr sync waves'],type='b',pch=20,col="blue")
# lines(centerinds,allres12[,'cor_>4yr sync NO3/waves'],type='b',pch=19,col="green")
# lines(centerinds,allres12[,'cor_>4yr sync NO3/waves'],type='b',lty="dashed",pch=20,col="blue",cex=.5)
# lines(centerinds,allres12[,'cor_>4yr sync NO3']+allres12[,'cor_>4yr sync waves']+allres12[,'cor_>4yr sync NO3/waves'],
#       type="b",col="yellow",pch=20)
# lines(centerinds,rep(0,length(centerinds)),type="l")
# put_loc_lines()
# dev.off()
# 
# #plot numbers of sites used - should be exactly the same as the corresponding kelp lag 4 plot
# pdf(paste0(resloc,"Numsitesused_v_centerind_kelplag12.pdf"))
# ylimits<-range(allres12[,'num locs'])
# plot(centerinds,allres12[,'num locs'],type='b',xlab="Center index",ylab="Number of sites used")
# lines(rep(PtConcInd,2),ylimits)
# lines(rep(SBInd,2),ylimits,lty="dashed")
# put_loc_lines()
# dev.off()

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

#So my three roughly comparable comparable regions are SBlocstouse, CC1locstouse, CC2locstouse

#***
#Now do the analysis for these three regions and save
#***

PlotForMSArgs<-list(cov_MainPlot_PanLabs=c("A","B","C"),
                    cov_MainPlot_short_ylim=NA, #means let the data determine the ylimits
                    cov_MainPlot_mid_ylim=NA,
                    cov_MainPlot_long_ylim=NA,
                    ExplainDirectMoran_fP1_short_ylims=c(0.0041,0.022),
                    ExplainDirectMoran_fP1_mid_ylims=c(0.00514,0.0136),
                    ExplainDirectMoran_fP1_long_ylims=c(0.00526,0.0126),
                    ExplainDirectMoran_Se1e1_short_ylims=c(0,95),
                    ExplainDirectMoran_Se1e1_mid_ylims=c(1.5,16),
                    ExplainDirectMoran_Se1e1_long_ylims=c(2,38),
                    ExplainDirectMoran_fP2_short_ylims=c(.01,.08),
                    ExplainDirectMoran_fP2_mid_ylims=c(.01,.08),
                    ExplainDirectMoran_fP2_long_ylims=c(.01,.08),
                    ExplainDirectMoran_Se2e2_short_ylims=c(0,8),
                    ExplainDirectMoran_Se2e2_mid_ylims=c(0.06,0.53),
                    ExplainDirectMoran_Se2e2_long_ylims=c(0.05,0.9),
                    ExplainDirectMoran_fB_short_ylims=c(0.4,3),
                    ExplainDirectMoran_fB_mid_ylims=c(0.8,1.65),
                    ExplainDirectMoran_fB_long_ylims=c(1.2,3),
                    ExplainDirectMoran_PanLabs=LETTERS[1:15],
                    ExplainInteractingMoran_fP1_short_ylims=c(0.064,0.15),
                    ExplainInteractingMoran_fP1_mid_ylims=c(0.0716,0.117),
                    ExplainInteractingMoran_fP1_long_ylims=c(0.0725,0.1125),
                    ExplainInteractingMoran_fP2_short_ylims=c(0.12,0.25),
                    ExplainInteractingMoran_fP2_mid_ylims=c(0.12,0.25),
                    ExplainInteractingMoran_fP2_long_ylims=c(0.12,0.25),
                    ExplainInteractingMoran_fP1ConjfP2_short_ylims=c(0.014,0.032),
                    ExplainInteractingMoran_fP1ConjfP2_mid_ylims=c(0.0139,0.0215),
                    ExplainInteractingMoran_fP1ConjfP2_long_ylims=c(0.01375,0.02),
                    ExplainInteractingMoran_Se1e2_short_ylims=c(0,30),
                    ExplainInteractingMoran_Se1e2_mid_ylims=c(0.04,1.4),
                    ExplainInteractingMoran_Se1e2_long_ylims=c(0.2,4),
                    ExplainInteractingMoran_Prod_short_ylims=c(0,0.7),
                    ExplainInteractingMoran_Prod_mid_ylims=c(0,0.028),
                    ExplainInteractingMoran_Prod_long_ylims=c(0.002,0.08),
                    ExplainInteractingMoran_PanLabs=LETTERS[1:15])
Args<-list(locstouse=CC1locstouse,lags=c(4,1,0),frg=c(0,.5),fnpre="RegionalAnalysis_CentCal1_KelpLag4",plottype="pdf",PlotForMSArgs=PlotForMSArgs)
res_CentCal1_KelpLag4<-do_analysis(Args)
saveRDS(res_CentCal1_KelpLag4,paste0(resloc,"res_CentCal1_KelpLag4.Rds"))

Args<-list(locstouse=CC1locstouse,lags=c(8,1,0),frg=c(0,.5),fnpre="RegionalAnalysis_CentCal1_KelpLag8",plottype="pdf")
res_CentCal1_KelpLag8<-do_analysis(Args)
saveRDS(res_CentCal1_KelpLag8,paste0(resloc,"res_CentCal1_KelpLag8.Rds"))

Args<-list(locstouse=CC1locstouse,lags=c(12,1,0),frg=c(0,.5),fnpre="RegionalAnalysis_CentCal1_KelpLag12",plottype="pdf")
res_CentCal1_KelpLag12<-do_analysis(Args)
saveRDS(res_CentCal1_KelpLag12,paste0(resloc,"res_CentCal1_KelpLag12.Rds"))

PlotForMSArgs<-list(cov_MainPlot_PanLabs=c("A","B","C"),
                    cov_MainPlot_short_ylim=NA, #means let the data determine the ylimits
                    cov_MainPlot_mid_ylim=NA,
                    cov_MainPlot_long_ylim=NA,
                    ExplainDirectMoran_fP1_short_ylims=c(0.0041,0.022),
                    ExplainDirectMoran_fP1_mid_ylims=c(0.00514,0.0136),
                    ExplainDirectMoran_fP1_long_ylims=c(0.00526,0.0126),
                    ExplainDirectMoran_Se1e1_short_ylims=c(0,95),
                    ExplainDirectMoran_Se1e1_mid_ylims=c(1.5,16),
                    ExplainDirectMoran_Se1e1_long_ylims=c(2,38),
                    ExplainDirectMoran_fP2_short_ylims=c(.01,.08),
                    ExplainDirectMoran_fP2_mid_ylims=c(.01,.08),
                    ExplainDirectMoran_fP2_long_ylims=c(.01,.08),
                    ExplainDirectMoran_Se2e2_short_ylims=c(0,8),
                    ExplainDirectMoran_Se2e2_mid_ylims=c(0.06,0.53),
                    ExplainDirectMoran_Se2e2_long_ylims=c(0.05,0.9),
                    ExplainDirectMoran_fB_short_ylims=c(0.4,3),
                    ExplainDirectMoran_fB_mid_ylims=c(0.8,1.65),
                    ExplainDirectMoran_fB_long_ylims=c(1.2,3),
                    ExplainDirectMoran_PanLabs=LETTERS[1:15],
                    ExplainInteractingMoran_fP1_short_ylims=c(0.064,0.15),
                    ExplainInteractingMoran_fP1_mid_ylims=c(0.0716,0.117),
                    ExplainInteractingMoran_fP1_long_ylims=c(0.0725,0.1125),
                    ExplainInteractingMoran_fP2_short_ylims=c(0.12,0.25),
                    ExplainInteractingMoran_fP2_mid_ylims=c(0.12,0.25),
                    ExplainInteractingMoran_fP2_long_ylims=c(0.12,0.25),
                    ExplainInteractingMoran_fP1ConjfP2_short_ylims=c(0.014,0.032),
                    ExplainInteractingMoran_fP1ConjfP2_mid_ylims=c(0.0139,0.0215),
                    ExplainInteractingMoran_fP1ConjfP2_long_ylims=c(0.01375,0.02),
                    ExplainInteractingMoran_Se1e2_short_ylims=c(0,30),
                    ExplainInteractingMoran_Se1e2_mid_ylims=c(0.04,1.4),
                    ExplainInteractingMoran_Se1e2_long_ylims=c(0.2,4),
                    ExplainInteractingMoran_Prod_short_ylims=c(0,0.7),
                    ExplainInteractingMoran_Prod_mid_ylims=c(0,0.028),
                    ExplainInteractingMoran_Prod_long_ylims=c(0.002,0.08),
                    ExplainInteractingMoran_PanLabs=LETTERS[1:15])
Args<-list(locstouse=CC2locstouse,lags=c(4,1,0),frg=c(0,.5),fnpre="RegionalAnalysis_CentCal2_KelpLag4",plottype="pdf",PlotForMSArgs=PlotForMSArgs)
res_CentCal2_KelpLag4<-do_analysis(Args)
saveRDS(res_CentCal2_KelpLag4,paste0(resloc,"res_CentCal2_KelpLag4.Rds"))

Args<-list(locstouse=CC2locstouse,lags=c(8,1,0),frg=c(0,.5),fnpre="RegionalAnalysis_CentCal2_KelpLag8",plottype="pdf")
res_CentCal2_KelpLag8<-do_analysis(Args)
saveRDS(res_CentCal2_KelpLag8,paste0(resloc,"res_CentCal2_KelpLag8.Rds"))

Args<-list(locstouse=CC2locstouse,lags=c(12,1,0),frg=c(0,.5),fnpre="RegionalAnalysis_CentCal2_KelpLag12",plottype="pdf")
res_CentCal2_KelpLag12<-do_analysis(Args)
saveRDS(res_CentCal2_KelpLag12,paste0(resloc,"res_CentCal2_KelpLag12.Rds"))

PlotForMSArgs<-list(cov_MainPlot_PanLabs=c("D","E","F"),
                    cov_MainPlot_short_ylim=NA, #means let the data determine the ylimits
                    cov_MainPlot_mid_ylim=NA,
                    cov_MainPlot_long_ylim=NA,
                    ExplainDirectMoran_fP1_short_ylims=c(0.0041,0.022),
                    ExplainDirectMoran_fP1_mid_ylims=c(0.00514,0.0136),
                    ExplainDirectMoran_fP1_long_ylims=c(0.00526,0.0126),
                    ExplainDirectMoran_Se1e1_short_ylims=c(0,95),
                    ExplainDirectMoran_Se1e1_mid_ylims=c(1.5,16),
                    ExplainDirectMoran_Se1e1_long_ylims=c(2,38),
                    ExplainDirectMoran_fP2_short_ylims=c(.01,.08),
                    ExplainDirectMoran_fP2_mid_ylims=c(.01,.08),
                    ExplainDirectMoran_fP2_long_ylims=c(.01,.08),
                    ExplainDirectMoran_Se2e2_short_ylims=c(0,8),
                    ExplainDirectMoran_Se2e2_mid_ylims=c(0.06,0.53),
                    ExplainDirectMoran_Se2e2_long_ylims=c(0.05,0.9),
                    ExplainDirectMoran_fB_short_ylims=c(0.4,3),
                    ExplainDirectMoran_fB_mid_ylims=c(0.8,1.65),
                    ExplainDirectMoran_fB_long_ylims=c(1.2,3),
                    ExplainDirectMoran_PanLabs=c(LETTERS[16:26],"AA","BB","CC","DD"),
                    ExplainInteractingMoran_fP1_short_ylims=c(0.064,0.15),
                    ExplainInteractingMoran_fP1_mid_ylims=c(0.0716,0.117),
                    ExplainInteractingMoran_fP1_long_ylims=c(0.0725,0.1125),
                    ExplainInteractingMoran_fP2_short_ylims=c(0.12,0.25),
                    ExplainInteractingMoran_fP2_mid_ylims=c(0.12,0.25),
                    ExplainInteractingMoran_fP2_long_ylims=c(0.12,0.25),
                    ExplainInteractingMoran_fP1ConjfP2_short_ylims=c(0.014,0.032),
                    ExplainInteractingMoran_fP1ConjfP2_mid_ylims=c(0.0139,0.0215),
                    ExplainInteractingMoran_fP1ConjfP2_long_ylims=c(0.01375,0.02),
                    ExplainInteractingMoran_Se1e2_short_ylims=c(0,30),
                    ExplainInteractingMoran_Se1e2_mid_ylims=c(0.04,1.4),
                    ExplainInteractingMoran_Se1e2_long_ylims=c(0.2,4),
                    ExplainInteractingMoran_Prod_short_ylims=c(0,0.7),
                    ExplainInteractingMoran_Prod_mid_ylims=c(0,0.028),
                    ExplainInteractingMoran_Prod_long_ylims=c(0.002,0.08),
                    ExplainInteractingMoran_PanLabs=c(LETTERS[16:26],"AA","BB","CC","DD"))
Args<-list(locstouse=SBlocstouse,lags=c(4,1,0),frg=c(0,.5),fnpre="RegionalAnalysis_SoCal_KelpLag4",plottype="pdf",PlotForMSArgs=PlotForMSArgs)
res_SoCal_KelpLag4<-do_analysis(Args)
saveRDS(res_SoCal_KelpLag4,paste0(resloc,"res_SoCal_KelpLag4.Rds"))

Args<-list(locstouse=SBlocstouse,lags=c(8,1,0),frg=c(0,.5),fnpre="RegionalAnalysis_SoCal_KelpLag8",plottype="pdf")
res_SoCal_KelpLag8<-do_analysis(Args)
saveRDS(res_SoCal_KelpLag8,paste0(resloc,"res_SoCal_KelpLag8.Rds"))

Args<-list(locstouse=SBlocstouse,lags=c(12,1,0),frg=c(0,.5),fnpre="RegionalAnalysis_SoCal_KelpLag12",plottype="pdf")
res_SoCal_KelpLag12<-do_analysis(Args)
saveRDS(res_SoCal_KelpLag12,paste0(resloc,"res_SoCal_KelpLag12.Rds"))

#DONE ***DAN: Done: Write some code that tracks and explains the contributions of the parts of the interaction
#term in the spectral equation. We want to see the phase relationship between waves and NO3, and the difference
#between the phase delays of the effects of these things on kelp, and how that produces different interaction
#effects at annual and >4yr timescales in central and southern CA. Figure out how to display this.

#Next steps: Find ways to picture the lag or the cross spectrum (or both) of the relationship between waves and 
#NO3 in the three regions. I already have some code somewhere that just looks at the average annual fluctuation 
#and shows a peak for waves in the winter and a peak for NO3 in the spring. Is this still true in So Cal? Are waves
#or nitrates less synchronous in So Cal? Or less influential on kelp? Or is just that the interaction effects have 
#disappeared because of shifting lags of effects? You can dig into So Cal more I think, start by looking closely
#again at the plots I already generated for So Cal. I will want to use elementary demonstrations of as many things
#as possible.

#Is the synchronizing effect of waves on kelp in SoCal less because waves are less synchronous, or less influential,
#or both? I can display the prefactor and the spectral term (real-parted and summed over off-diagonal entries, probably)
#for NO3 and for waves, and compare those also acrross the three regions. Likewise we can look at fP1*Conj(fP2) 
#and also Se1e2 (appropriately summed across non-diagonal entries) and compare them across the three regions to
#get a sense whether we are talking about decreased cross-location synchrony between NO3 and waves (relative to the
#other regions) or modified influence of these factors on kelp.

#DECIDED NO Should my various component plots be using the prefactor 1/|fB|^2?

#Thoughts for next time I get back to work on this:
#1) Small presentational thing: Since you are not partitioning the timescale range anyway (you don’t really do anything with 2-4 year 
#timescales), and since the purpose of this paper is to illustrate concepts about interacting Moran effects and
#not to explain everything about kelp, might as well focus on the immediate neighborhood of the annual timescale 
#(where you can see the peaks) and the longest timescale (greater than 8 or 16 years, no specific need to start at 4). 
#Results are cleanest there. Put grey rectangles on plots covering those areas. You can still use the dividers
#at 2 and 4 years for the panel separation, no problem, since those are just to allow different y-axis extents.
#Though there may be an argument for combining the 2-4 and >4 yr range panels into one panel, since y-axis extents 
#tend to be similar for those. I DECIDED TO FOCUS ON THOSE RANGES, AS INDICATED, BUT NOT TO PUT THE RECTANGLES.

