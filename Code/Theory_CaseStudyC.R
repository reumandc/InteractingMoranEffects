#This script implements theoretical case study C

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): matrixcalc, mvtnorm
source("SpectralTools.R")
source("WhichBreak.R")
source("plotColLine.R")

#***
#Location for storing results and other prep
#***

resloc<-"../Results/Theory/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***
#The main function for doing both case studies
#***

#***
#The main function for doing both case studies
#***

#This function calculates the analytic results and plots them for case study C.
#The function is actually more general than the planned case study that will be 
#displayed in the paper.
#
#Args
#b            The parameters (b_1,...,b_n) of the model
#p1           The parameters (p_0^{(1)},...,p_{m_1}^{(1)}) of the model
#p2           The parameters (p_0^{(2)},...,p_{m_1}^{(2)}) of the model
#N            The number of habitat patches
#ns           Parameters for the noise. A length-7 vector with some constraints.
#               See the code or the paper for information on how the noise is
#               constructed based on these parameters. Using notation from the paper,
#               the entries are c1, c2, l1, l2, w11, w22, w12.
#PlotName     A file name of the saved plots, including path, excluding file extension.
#PanLabs      A vector of panel labels (typically c("A","B","C","D"), or similar)
#PanLabs2     Another similar, this one just 3 long, for the sup mat figure showing the 
#               rhos
#SyncYAx      A length 2 vector with the y axis extent for the synchrony plot. This is 
#               just for the publication-quality plotting, set to NA if you want to 
#               have this determined automatically.
#RelSyncYAx   Similar but for the relative synchrony plot.
#XAxOn        Default TRUE means plot the x axis label and numbers, FALSE means don't,
#               with the size of the plot correspondingly adjusted.
#
#Output - none, but plots saved.
#
#Notes
#Very little error checking, so watch out!
#
TheoryCaseC<-function(b,p1,p2,N,ns,PlotName,PanLabs,PanLabs2,SyncYAx=NA,RelSyncYAx=NA,XAxOn=TRUE)
{
  #construct Sig, a covariance matrix used in construction of the noise, basically 
  #just to find out if it is positive definite
  w11<-ns[5]
  Sig11<-matrix(w11,N,N)
  diag(Sig11)<-1
  w22<-ns[6]
  Sig22<-matrix(w22,N,N)
  diag(Sig22)<-1
  w12<-ns[7]
  Sig12<-matrix(w12,N,N)
  SigL<-rbind(Sig11,Sig12)
  SigR<-rbind(Sig12,Sig22)
  Sig<-cbind(SigL,SigR)
  if (!matrixcalc::is.positive.definite(Sig))
  {
    stop("Error in TheoryCaseC: bad values for ns, resulted in non-positive definite covariance matrix")
  }
  
  #construct fB, fP1, fP2
  Tx<-256 #These are analyic results, so we don't have a time series length,
  #but we're plotting as though we had time series of this length. Make this 
  #a power of 2. Plotting not guaranteed to work corectly if you change this. 
  freq<-freq<-(1:(Tx-1))/Tx
  freq<-freq[freq<=0.5]
  mu<-exp(-2*pi*complex(real=0,imaginary=1)*freq)
  fP1<-complex(real=numeric(length(freq)),imaginary=numeric(length(freq)))
  for (counter in 1:length(p1))
  {
    fP1<-fP1+p1[counter]*mu^(counter-1)
  }
  fP2<-complex(real=numeric(length(freq)),imaginary=numeric(length(freq)))
  for (counter in 1:length(p2))
  {
    fP2<-fP2+p2[counter]*mu^(counter-1)
  }
  fB<-complex(real=numeric(length(freq))+1,imaginary=numeric(length(freq)))
  for (counter in 1:length(b))
  {
    fB<-fB-b[counter]*mu^counter
  }
  
  #now construct the rhos - simple, based on the math in the sup mat of the paper
  c1<-ns[1]
  c2<-ns[2]
  l1<-ns[3]
  l2<-ns[4]
  rho11<-w11*(Mod(1/(1-c1*mu-c2*mu^2)))^2
  rho22<-w22*(Mod(1/(1-c1*mu-c2*mu^2)))^2
  iota<-complex(real=0,imaginary=1)
  rho12<-w12*(Mod(1/(1-c1*mu-c2*mu^2)))^2*exp(-2*pi*iota*freq*(l1-l2))
  
  #construct the quantities to be plotted - starting with the direct Moran effect of epsilon^(1)
  direct1<-rho11*((Mod(fP1))^2)/((Mod(fB))^2)
  direct1_wopd<-rho11*((Mod(fP1))^2)
  
  #construct the quantities to be plotted - now the direct Moran effect of epsilon^(2)
  direct2<-rho22*((Mod(fP2))^2)/((Mod(fB))^2)
  direct2_wopd<-rho22*((Mod(fP2))^2)
  
  #construct the quantities to be plotted - now the interaction term
  interact<-2*Re(fP1*Conj(fP2)*rho12)/((Mod(fB))^2)
  interact_wopd<-2*Re(fP1*Conj(fP2)*rho12)
  
  #Now the terms inside the Re(.) operator - first fP1*conj(fP2), then rho12
  interact_p1<-fP1*Conj(fP2)
  interact_p2<-rho12
  
  #***Generate the main figure you need for this case study
  
  #dimensions for the figure to be generated, inches
  totwd<-6.5
  xaxht<-0.45
  yaxwd<-0.45
  gap<-0.1
  panwd<-(totwd-2*gap-3*yaxwd)/4
  panht<-panwd
  if (XAxOn)
  {
    totht<-xaxht+panht+gap
  } else
  {
    totht<-gap+panht+gap
  }
  lnwid<-2.5
  cexmath<-0.9
  pdf(file=paste0(PlotName,".pdf"),width=totwd,height=totht)
  
  #prep for plotting
  timescales<-1/freq
  l2timescales<-log2(timescales)
  
  #technical stuff for the circular colorbar for phases fo the below panels
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
  
  #prep for panel 1  
  if (XAxOn)
  {
    par(fig=c((yaxwd)/totwd,
              (yaxwd+panwd)/totwd,
              (xaxht)/totht,
              (xaxht+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
  } else
  {
    par(fig=c((yaxwd)/totwd,
              (yaxwd+panwd)/totwd,
              (gap)/totht,
              (gap+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
  }
  
  #Now do the panel that shows the direct Moran contributions to
  #synchrony and the interaction effect
  xlimits<-range(l2timescales)
  if (identical(SyncYAx,NA))
  {
    ylimits<-range(direct1,direct2,interact)
    ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  } else
  {
    ylimits<-SyncYAx  
  }
  
  plot(l2timescales,direct1,type="l",ylim=ylimits,col="green",xaxt="n",ylab="",xlab="",xaxs="i",lwd=lnwid)
  lines(l2timescales,direct2,type="l",col="blue",lwd=lnwid)
  lines(l2timescales,interact,type="l",col="green",lwd=lnwid)
  lines(l2timescales,interact,type="l",lty="dashed",col="blue",lwd=lnwid)
  
  mtext("Synch. contrib.",2,1.2)
  l2tslocs<-seq(from=min(l2timescales),by=1,to=max(l2timescales))
  if (XAxOn)
  {
    graphics::axis(1,at=l2tslocs,labels=c("0.5","","2","","8","","32",""))
    mtext("Timescale, years",1,1.2)
  } else
  {
    graphics::axis(1,at=l2tslocs,labels=FALSE)
  }    
  lines(range(l2timescales),c(0,0),type="l",lty="dotted")
  text(xlimits[1],ylimits[2],PanLabs[1],cex=1.2,adj=c(-0.1,1))
  
  #prep for panel 2
  if (XAxOn)
  {
    par(fig=c((2*yaxwd+panwd)/totwd,
              (2*yaxwd+2*panwd)/totwd,
              (xaxht)/totht,
              (xaxht+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  } else
  {
    par(fig=c((2*yaxwd+panwd)/totwd,
              (2*yaxwd+2*panwd)/totwd,
              (gap)/totht,
              (gap+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  }
  
  #Make the panel that shows the phases and magnitude of the 
  #some terms inside the Re(.) for the interaction effects. Looking at fP1, conj(fP2),
  #and fP2*conj(fP2) here.
  ylimits<-range(Mod(fP1),Mod(fP2),Mod(fP1*fP2),Mod(rho12),Mod(fP1*fP2*rho12))
  di<-diff(ylimits)
  ylimits[2]<-ylimits[2]+.5*di
  ylimits[1]<-ylimits[1]-.2*di
  
  y<-Mod(fP1)
  h<-Arg(fP1)
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales,y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xaxs="i",xaxt="n")
  text(l2timescales[length(l2timescales)-10],y[length(l2timescales)-10],latex2exp::TeX("$f_{P^{(1)}}$"),
       cex=cexmath,adj=c(0,0))
  
  y<-Mod(fP2)
  h<-Arg(fP2)
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  linesColLine(x=l2timescales,y=y,cols=cols,lnwid=lnwid)
  text(l2timescales[length(l2timescales)-10],y[length(l2timescales)-10],latex2exp::TeX("$f_{P^{(2)}}$"),
       cex=cexmath,adj=c(0,1))
  
  y<-Mod(fP1*Conj(fP2))
  h<-Arg(fP1*Conj(fP2))
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  linesColLine(x=l2timescales,y=y,cols=cols,lnwid=lnwid)  
  text(l2timescales[length(l2timescales)-10],y[length(l2timescales)-10],
       latex2exp::TeX("$f_{P^{(1)}} \\bar{f_{P^{(2)}}}$"),
       cex=cexmath,adj=c(0,0))
  
  mtext("Magnitude",2,1.2)
  l2tslocs<-seq(from=min(l2timescales),by=1,to=max(l2timescales))
  if (XAxOn)
  {
    graphics::axis(1,at=l2tslocs,labels=c("0.5","","2","","8","","32",""))
    mtext("Timescale, years",1,1.2)
  } else
  {
    graphics::axis(1,at=l2tslocs,labels=FALSE)
  }    
  text(xlimits[1],ylimits[2],PanLabs[2],cex=1.2,adj=c(-0.1,1))
  
  #make the colorbar
  bds<-seq(from=.6*xlimits[1]+.4*xlimits[2],to=.1*xlimits[1]+.9*xlimits[2],
           length.out=length(breaks))
  rect(xleft=bds[1:(length(bds)-1)],
       ybottom=.1*ylimits[1]+.9*ylimits[2],
       xright=bds[2:length(bds)],
       ytop=ylimits[2],
       col=colbar,
       border=NA)
  ytags<-c(.1*ylimits[1]+.9*ylimits[2],.15*ylimits[1]+.85*ylimits[2],.175*ylimits[1]+.825*ylimits[2])
  lines(rep(.6*xlimits[1]+.4*xlimits[2],2),ytags[1:2])
  text(.6*xlimits[1]+.4*xlimits[2],ytags[3],latex2exp::TeX("$-\\pi$"),adj=c(.5,1))
  lines(rep((.6-.25)*xlimits[1]+(.4+.25)*xlimits[2],2),ytags[1:2])
  text((.6-.25)*xlimits[1]+(.4+.25)*xlimits[2],ytags[3],"0",adj=c(.5,1))
  lines(rep(.1*xlimits[1]+.9*xlimits[2],2),ytags[1:2])
  text(.1*xlimits[1]+.9*xlimits[2],ytags[3],latex2exp::TeX("$\\pi$"),adj=c(.5,1))
  
  #prep for panel 3
  if (XAxOn)
  {
    par(fig=c((2*yaxwd+2*panwd+gap)/totwd,
              (2*yaxwd+3*panwd+gap)/totwd,
              (xaxht)/totht,
              (xaxht+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  } else
  {
    par(fig=c((2*yaxwd+2*panwd+gap)/totwd,
              (2*yaxwd+3*panwd+gap)/totwd,
              (gap)/totht,
              (gap+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  }
  
  #make the panel that continues to show how the terms inside the
  #Re(.) combine. 
  y<-Mod(fP1*Conj(fP2))
  h<-Arg(fP1*Conj(fP2))
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales,y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xaxs="i",xaxt="n",yaxt="n")
  text(l2timescales[length(l2timescales)-10],y[length(l2timescales)-10],latex2exp::TeX("$f_{P^{(1)}} \\bar{f_{P^{(2)}}}$"),
       cex=cexmath,adj=c(0,0))
  
  y<-rep(Mod(rho12),length(l2timescales))
  h<-rep(Arg(rho12),length(l2timescales))
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  linesColLine(x=l2timescales,y=y,cols=cols,lnwid=lnwid)
  text(l2timescales[length(l2timescales)-10],y[length(l2timescales)-10],
       latex2exp::TeX("$\\rho_{\\epsilon^{(1)} \\epsilon^{(2)}}$"),
       cex=cexmath,adj=c(0,1.2))
  
  y<-Mod(fP1*Conj(fP2)*rho12)
  h<-Arg(fP1*Conj(fP2)*rho12)
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  linesColLine(x=l2timescales,y=y,cols=cols,lnwid=lnwid)
  text(l2timescales[length(l2timescales)-90],y[length(l2timescales)-90],
       latex2exp::TeX("$f_{P^{(1)}} \\bar{f_{P^{(2)}}} \\rho_{\\epsilon^{(1)} \\epsilon^{(2)}}$"),
       cex=cexmath,adj=c(0,0))
  
  if (XAxOn)
  {
    graphics::axis(1,at=l2tslocs,labels=c("0.5","","2","","8","","32",""))
    mtext("Timescale, years",1,1.2)
  } else
  {
    graphics::axis(1,at=l2tslocs,labels=FALSE)
  }    
  graphics::axis(2,labels=FALSE)
  text(xlimits[1],ylimits[2],PanLabs[3],cex=1.2,adj=c(-0.1,1))
  
  #prep for panel 4
  if (XAxOn)
  {
    par(fig=c((3*yaxwd+3*panwd+gap)/totwd,
              (3*yaxwd+4*panwd+gap)/totwd,
              (xaxht)/totht,
              (xaxht+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  } else
  {
    par(fig=c((3*yaxwd+3*panwd+gap)/totwd,
              (3*yaxwd+4*panwd+gap)/totwd,
              (gap)/totht,
              (gap+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  }
  
  #make the panel that shows the direct Moran contributions to
  #synchrony and the interaction effect, but without the modifications
  #made by the fB term
  
  if (identical(RelSyncYAx,NA))
  {
    ylimits<-range(direct1_wopd,direct2_wopd,interact_wopd)
    ylimits[2]<-ylimits[2]+.2*diff(ylimits)
  } else
  {
    ylimits<-RelSyncYAx  
  }
  
  plot(l2timescales,direct1_wopd,type="l",ylim=ylimits,col="green",xaxt="n",ylab="",xlab="",xaxs="i",lwd=lnwid)
  lines(l2timescales,direct2_wopd,type="l",col="blue",lwd=lnwid)
  lines(l2timescales,interact_wopd,type="l",col="green",lwd=lnwid)
  lines(l2timescales,interact_wopd,type="l",lty="dashed",col="blue",lwd=lnwid)
  mtext("Rel. synch. contrib.",2,1.2)
  if (XAxOn)
  {
    graphics::axis(1,at=l2tslocs,labels=c("0.5","","2","","8","","32",""))
    mtext("Timescale, years",1,1.2)
  } else
  {
    graphics::axis(1,at=l2tslocs,labels=FALSE)
  }    
  lines(range(l2timescales),c(0,0),type="l",lty="dotted")
  text(xlimits[1],ylimits[2],PanLabs[4],cex=1.2,adj=c(-0.1,1))
  
  dev.off()
  
  #***For case study C, we also need a separate figure which just displace the rhos
  
  #dimensions for the figure to be generated, inches
  totwd<-6.5
  xaxht<-0.45
  yaxwd<-0.45
  gap<-0.1
  panwd<-(totwd-gap-3*yaxwd)/3
  panht<-panwd
  if (XAxOn)
  {
    totht<-xaxht+panht+gap
  } else
  {
    totht<-gap+panht+gap
  }
  lnwid<-2.5
  pdf(file=paste0(PlotName,"_rhos.pdf"),width=totwd,height=totht)
  
  ylimits<-range(Mod(rho11),Mod(rho22),Mod(rho12))
  ylimits[2]<-ylimits[2]+.5*diff(ylimits)
  
  #prep for panel 1 - rho11  
  if (XAxOn)
  {
    par(fig=c((yaxwd)/totwd,
              (yaxwd+panwd)/totwd,
              (xaxht)/totht,
              (xaxht+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
  } else
  {
    par(fig=c((yaxwd)/totwd,
              (yaxwd+panwd)/totwd,
              (gap)/totht,
              (gap+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
  }
  
  y<-Mod(rho11)
  h<-Arg(rho11)
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales,y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xaxs="i",xaxt="n",yaxt="n")

  if (XAxOn)
  {
    graphics::axis(1,at=l2tslocs,labels=c("0.5","","2","","8","","32",""))
    mtext("Timescale, years",1,1.2)
  } else
  {
    graphics::axis(1,at=l2tslocs,labels=FALSE)
  }    
  mtext(latex2exp::TeX("$\\rho_{\\epsilon^{(1)} \\epsilon^{(1)}}$"),2,1.2)
  graphics::axis(2,labels=TRUE)
  text(xlimits[1],ylimits[2],PanLabs2[1],cex=1.2,adj=c(-0.1,1))
  
  #make the colorbar
  bds<-seq(from=.6*xlimits[1]+.4*xlimits[2],to=.1*xlimits[1]+.9*xlimits[2],
           length.out=length(breaks))
  rect(xleft=bds[1:(length(bds)-1)],
       ybottom=.1*ylimits[1]+.9*ylimits[2],
       xright=bds[2:length(bds)],
       ytop=ylimits[2],
       col=colbar,
       border=NA)
  ytags<-c(.1*ylimits[1]+.9*ylimits[2],.15*ylimits[1]+.85*ylimits[2],.175*ylimits[1]+.825*ylimits[2])
  lines(rep(.6*xlimits[1]+.4*xlimits[2],2),ytags[1:2])
  text(.6*xlimits[1]+.4*xlimits[2],ytags[3],latex2exp::TeX("$-\\pi$"),adj=c(.5,1))
  lines(rep((.6-.25)*xlimits[1]+(.4+.25)*xlimits[2],2),ytags[1:2])
  text((.6-.25)*xlimits[1]+(.4+.25)*xlimits[2],ytags[3],"0",adj=c(.5,1))
  lines(rep(.1*xlimits[1]+.9*xlimits[2],2),ytags[1:2])
  text(.1*xlimits[1]+.9*xlimits[2],ytags[3],latex2exp::TeX("$\\pi$"),adj=c(.5,1))
  
  #prep for panel 2 - rho22
  if (XAxOn)
  {
    par(fig=c((2*yaxwd+panwd)/totwd,
              (2*yaxwd+2*panwd)/totwd,
              (xaxht)/totht,
              (xaxht+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  } else
  {
    par(fig=c((2*yaxwd+panwd)/totwd,
              (2*yaxwd+2*panwd)/totwd,
              (gap)/totht,
              (gap+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  }
  
  y<-Mod(rho22)
  h<-Arg(rho22)
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales,y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xaxs="i",xaxt="n",yaxt="n")
  
  if (XAxOn)
  {
    graphics::axis(1,at=l2tslocs,labels=c("0.5","","2","","8","","32",""))
    mtext("Timescale, years",1,1.2)
  } else
  {
    graphics::axis(1,at=l2tslocs,labels=FALSE)
  }    
  mtext(latex2exp::TeX("$\\rho_{\\epsilon^{(2)} \\epsilon^{(2)}}$"),2,1.2)
  graphics::axis(2,labels=TRUE)
  text(xlimits[1],ylimits[2],PanLabs2[2],cex=1.2,adj=c(-0.1,1))
  
  #prep for panel 3  
  if (XAxOn)
  {
    par(fig=c((3*yaxwd+2*panwd)/totwd,
              (3*yaxwd+3*panwd)/totwd,
              (xaxht)/totht,
              (xaxht+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  } else
  {
    par(fig=c((3*yaxwd+2*panwd)/totwd,
              (3*yaxwd+3*panwd)/totwd,
              (gap)/totht,
              (gap+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  }
  
  y<-Mod(rho12)
  h<-Arg(rho12)
  h[h==-pi]<-pi
  cols<-colbar[WhichBreak(breaks=breaks,h)]
  plotColLine(x=l2timescales,y=y,cols=cols,lnwid=lnwid,ylim=ylimits,xaxs="i",xaxt="n",yaxt="n")
  
  if (XAxOn)
  {
    graphics::axis(1,at=l2tslocs,labels=c("0.5","","2","","8","","32",""))
    mtext("Timescale, years",1,1.2)
  } else
  {
    graphics::axis(1,at=l2tslocs,labels=FALSE)
  }    
  mtext(latex2exp::TeX("$\\rho_{\\epsilon^{(1)} \\epsilon^{(2)}}$"),2,1.2)
  graphics::axis(2,labels=TRUE)
  text(xlimits[1],ylimits[2],PanLabs2[3],cex=1.2,adj=c(-0.1,1))
  
  dev.off()
}

#Case C, parameter set 1
N<-5
ns<-c(0,-4/9,0,1,.6,.6,.3) #c1, c2, l1, l2, w11, w22, w12.
p1<-c(0,3)
p2<-c(1.5)
b<-c(.4)
PlotName<-paste0(resloc,"TheoryFigCaseC_v01")
PanLabs<-c("i","j","k","l")
PanLabs2<-c("a","b","c")
TheoryCaseC(b,p1,p2,N,ns,PlotName,PanLabs,PanLabs2,SyncYAx=c(-10,20),RelSyncYAx=c(-10,25),XAxOn=FALSE)

#Case study C, parameter set 2
ns<-c(0,-4/9,1,0,.6,.6,.3) #c1, c2, l1, l2, w11, w22, w12.
p1<-c(0,3)
p2<-c(1.5)
b<-c(.4)
PlotName<-paste0(resloc,"TheoryFigCaseC_v02")
PanLabs<-c("m","n","o","p")
PanLabs2<-c("d","e","f")
TheoryCaseC(b,p1,p2,N,ns,PlotName,PanLabs,PanLabs2,SyncYAx=c(-10,20),RelSyncYAx=c(-10,25),XAxOn=TRUE)
