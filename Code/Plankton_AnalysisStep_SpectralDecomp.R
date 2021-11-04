#Spectral analysis steps applied to plankton.

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

resloc<-"../Results/Plankton_SpectralDecompResults/"
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
#The main function for doing the analysis 
#***

#Uses the spectral mathematics to get information about fractions of synchrony in pci explained by temperature, cal fin, 
#and interactions. Assumes all the variables created above are in the workspace when run.
#
#Args - a list with these named elements
#locstouse        Indices of locations to use. These refer to rows of the pci, temp and calfin data matrices above
#lags             Max regression lags to use, pci then temp then calfin. In years. So if you want the model to
#                   allow a pci lag effect of 1 year and temp and calfin lag effects from the current year 
#                   and last year, you input lags=c(1,1,1). Lags used will be 1 up to the value input for pci,
#                   and 0 up to the value input for temp, and the same for calfin.
#frg              The frequency range to use for display
#fnpre            A prefix for file names of pdf plots exported by the function, with path but without the
#                   extension "pdf" or "jpg"
#plottype         Use "jpg" to make jpg plots, "pdf" to make pdf plots
#PlotForMSArgs    Formatting arguments used for generating final, presentation-quality plots. If this is missing it 
#                   means don't generate those plots.
#
do_analysis<-function(Args)
{
  #**argument unpacking and error catching - global vars accessed here (pci)
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
  if (any(locstouse<1) || any(locstouse>dim(pci)[1]))
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
  
  #**filter the data to only the specified locations - global vars accessed here (pci, temp, calfin)
  pci_l<-pci[locstouse,]
  temp_l<-temp[locstouse,]
  calfin_l<-calfin[locstouse,]
  numlocs_l<-length(locstouse)

  #**do regression to get coefficients for fB, fP1, fP2
  
  #construct the regression formula object
  form_rhs_terms<-c(paste0("pci_l",1:lags[1]),paste0("temp_l",0:lags[2]),paste0("calfin_l",0:lags[3]))
  form<-as.formula(paste0("pci_l0~",paste0(form_rhs_terms,collapse="+"),"-1",sep=""))
  
  #construct the regression data frame
  maxlag<-max(lags)
  regdat<-as.data.frame(matrix(NA,(dim(pci_l)[1])*(dim(pci_l)[2]-maxlag),length(form_rhs_terms)+1))
  names(regdat)<-c("pci_l0",form_rhs_terms)
  for (counter in 0:lags[1]) #counter is the pci lag here
  {
    regdat[,counter+1]<-as.vector(pci_l[,(maxlag+1-counter):(lents-counter)])
  }
  for (counter in 0:lags[2]) #counter is the temp lag here
  {
    regdat[,counter+lags[1]+2]<-as.vector(temp_l[,(maxlag+1-counter):(lents-counter)])
  }
  for (counter in 0:lags[3]) #counter is the calfin lag here
  {
    regdat[,counter+lags[1]+1+lags[2]+1+1]<-as.vector(calfin_l[,(maxlag+1-counter):(lents-counter)])
  }
  
  #do the regression, extract coefficients
  mod<-lm(formula=form,data=regdat)
  cm<-coef(mod)
  
  #**get residuals, cut data to common length, and compute the necessary spectral matrices
  
  #cut down the data - global vars accessed here (lents)
  delta<-matrix(unname(residuals(mod)),numlocs_l,lents-maxlag)
  pci_l<-pci_l[,(maxlag+1):lents]
  temp_l<-temp_l[,(maxlag+1):lents]
  calfin_l<-calfin_l[,(maxlag+1):lents]
  
  #get spectral matrices - global vars accessed here (BiasVariance)
  #e1 is temp, e2 is calfin, e is combined temp and then calfin, w is pci
  Sww<-myspecmatbrill(pci_l,detrend=FALSE,cutsym=TRUE,forvar=FALSE,BiasVariance=BiasVariance)
  freq<-Sww$freq
  Sww<-Sww$spec
  S<-myspecmatbrill(rbind(temp_l,calfin_l,delta),detrend=FALSE,cutsym=TRUE,forvar=FALSE,BiasVariance=BiasVariance)
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
    fP1<-fP1+unname(cm[lags[1]+1+counter])*mu^counter
  }
  fP2<-complex(real=numeric(length(freq)),imaginary=numeric(length(freq)))
  for (counter in 0:lags[3])
  {
    fP2<-fP2+unname(cm[lags[1]+lags[2]+1+1+counter])*mu^counter
  }
  fB<-complex(real=numeric(length(freq))+1,imaginary=numeric(length(freq)))
  for (counter in 1:lags[1])
  {
    fB<-fB-unname(cm[counter])*mu^counter
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
  
  #**make plots - this is the "_MainPlot" plot that shows all the components of synchrony
  make_plots(freq,frg,totsync_cov,T1_avg_cov,T2_avg_cov,T3_avg_cov,T4_avg_cov,T5_avg_cov,T6_avg_cov,resloc,
             fnpre=paste0(fnpre,"_cov"),plottype)
  
  
  return(NA)
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
  
  #we are going to plot against log2 timescale, on two panels for different timescale ranges, so
  #prep for that
  timescales<-1/freq
  l2timescales<-log2(timescales)
  xlimits<-range(l2timescales)
  xlimits_short<-c(xlimits[1],log2(4))
  xlimits_long<-c(log2(4),xlimits[2])
  
  #panel layout quantities
  gap<-.15
  yaxnumwd<-.4
  yaxlabwd<-.2
  yaxwd<-yaxnumwd+yaxlabwd
  xaxht<-yaxwd
  panht<-3
  panwd_short<-2
  panwd_long<-2
  totpanwd<-panwd_short+panwd_long
  totht<-xaxht+panht+gap
  totwd<-yaxwd+1*yaxnumwd+totpanwd+gap
  
  #get the plotting device
  if (plottype=="jpg")
  {
    jpeg(paste0(resloc,fnpre,"_MainPlot.jpg"),quality=95,width=totwd,height=totht,units="in",res=300)
  }
  if (plottype=="pdf")
  {
    pdf(paste0(resloc,fnpre,"_MainPlot.pdf"),width=totwd,height=totht)
  }
  
  #the panel for short timescales
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
  graphics::axis(1,at=log2(c(2,4)),labels=c(2,4))
  lines(l2timescales[inds],T1_avg[inds],type="l",col="green") #direct effects of temp
  lines(l2timescales[inds],T2_avg[inds],type="l",col="blue") #direct effects of calfin
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="solid",col="green") 
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="dashed",col="blue") #interactions between temp and calfin
  lines(l2timescales[inds],sumexpl[inds],type="l",col="red",lty="solid") #explained stuff - sum of all the terms relating only to temp and calfin (direct temp and calfin effects, and interactions between those two)
  lines(l2timescales[inds],sumterms[inds],type="l",lty="dashed") #total of all terms, should approx equal the plot of totsync
  #lines(l2timescales[inds],sumunexpl[inds],type="l",lty="solid",col="red") #unexplained
  lines(l2timescales[inds],rep(0,length(timescales[inds])))
  mtext("Timescale, years",1,1.7)
  
  #the panel for long (>4 yr) timescales
  par(fig=c((yaxwd+panwd_short+yaxnumwd)/totwd,
            (yaxwd+panwd_short+yaxnumwd+panwd_long)/totwd,
            (xaxht)/totht,
            (xaxht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,0.75,0),new=TRUE)
  inds<-which(l2timescales>=xlimits_long[1] & l2timescales<=xlimits_long[2])
  ylimits<-range(totsync[inds],T1_avg[inds],T2_avg[inds],T4_avg[inds],sumexpl[inds],sumterms[inds])
  plot(l2timescales[inds],totsync[inds],type="l",xaxt="n",col="black",
       ylim=ylimits,xlim=xlimits_long,xaxs="i")
  graphics::axis(1,at=log2(c(4,8,16,32)),labels=c(4,8,16,32))
  lines(l2timescales[inds],T1_avg[inds],type="l",col="green") #direct effects of temp
  lines(l2timescales[inds],T2_avg[inds],type="l",col="blue") #direct effects of calfin
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="solid",col="green") 
  lines(l2timescales[inds],T4_avg[inds],type="l",lty="dashed",col="blue") #interactions between temp and calfin
  lines(l2timescales[inds],sumexpl[inds],type="l",col="red",lty="solid") #explained stuff 
  lines(l2timescales[inds],sumterms[inds],type="l",lty="dashed") #total of all terms, should approx equal the plot of totsync
  #lines(l2timescales[inds],sumunexpl[inds],type="l",lty="solid",col="red") #unexplained
  lines(l2timescales[inds],rep(0,length(timescales[inds])))
  mtext("Timescale, years",1,1.7)
  
  dev.off()
  
  return(NULL)
}

#***
#Now do the analysis and save
#***

Args<-list(locstouse=1:(dim(pci)[1]),lags=c(1,1,0),frg=c(0,.5),fnpre="AllLocs_Lags1_1_0",plottype="pdf",PlotForMSArgs=NA)
res_AllLocs_Lags1_1_0<-do_analysis(Args)

Args<-list(locstouse=1:(dim(pci)[1]),lags=c(1,1,1),frg=c(0,.5),fnpre="AllLocs_Lags1_1_1",plottype="pdf",PlotForMSArgs=NA)
res_AllLocs_Lags1_1_1<-do_analysis(Args)

Args<-list(locstouse=1:(dim(pci)[1]),lags=c(2,2,2),frg=c(0,.5),fnpre="AllLocs_Lags2_2_2",plottype="pdf",PlotForMSArgs=NA)
res_AllLocs_Lags2_2_2<-do_analysis(Args)

Args<-list(locstouse=1:(dim(pci)[1]),lags=c(5,5,5),frg=c(0,.5),fnpre="AllLocs_Lags5_5_5",plottype="pdf",PlotForMSArgs=NA)
res_AllLocs_Lags5_5_5<-do_analysis(Args)

Args<-list(locstouse=1:(dim(pci)[1]),lags=c(10,10,10),frg=c(0,.5),fnpre="AllLocs_Lags10_10_10",plottype="pdf",PlotForMSArgs=NA)
res_AllLocs_Lags10_10_10<-do_analysis(Args)

Args<-list(locstouse=1:(dim(pci)[1]),lags=c(5,1,0),frg=c(0,.5),fnpre="AllLocs_Lags5_1_0",plottype="pdf",PlotForMSArgs=NA)
res_AllLocs_Lags5_1_0<-do_analysis(Args)


#saveRDS(res_CentCal1_KelpLag4,paste0(resloc,"res_CentCal1_KelpLag4.Rds"))

