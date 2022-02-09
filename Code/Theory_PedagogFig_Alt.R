#This script makes an alternative pedagogical figure, may turn out to be an improvement. No modelling here, 
#just presentation. For the early part of the paper.

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): shape, latex2exp

#***
#Location for storing results and other prep
#***

resloc<-"../Results/Theory/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

set.seed(102)

#***
#Make the figure
#***

#dimensions for the figure to be generated, inches
titlespace<-.25
totwd<-3.5
gap<-0.25
xaxht<-0.45
yaxwd<-0.25
panwd<-(totwd-gap-yaxwd)
panht<-panwd/2.1
totht<-(xaxht+panht+gap+panht+titlespace)*2

pdf(file=paste0(resloc,"PedagogFig_Alt.pdf"),height=totht,width=totwd)

#***part 1, which shows synergistic effects

#get all the quantities to plot
tim<-seq(from=0,to=80,length.out=1001)
period<-20

ln<-3
phase1<-(-1.5)
phase2<-phase1+2*pi*ln/period

phase_jitters_1<-2*pi*(2*runif(3)-1)/period
phase_jitters_2<-2*pi*(3*runif(3)-1)/period

vertoffset<-8

envvar_1<-matrix(NA,3,length(tim))
envvar_2<-matrix(NA,3,length(tim))
for (counter in 1:(dim(envvar_1)[1]))
{
  envvar_1[counter,]<-sin(2*pi*tim/period+phase1+phase_jitters_1[counter])+(counter-1)*vertoffset
  envvar_2[counter,]<-sin(2*pi*tim/period+phase2+phase_jitters_2[counter])+(counter-1)*vertoffset
}

sinpks_1_x<-matrix(NA,3,4)
sinpks_1_x[,1]<-(pi/2-phase1-phase_jitters_1)*period/2/pi
sinpks_2_x<-matrix(NA,3,4)
sinpks_2_x[,1]<-(pi/2-phase2-phase_jitters_2)*period/2/pi
for (counter in 2:4)
{
  sinpks_1_x[,counter]<-sinpks_1_x[,1]+period*(counter-1)
  sinpks_2_x[,counter]<-sinpks_2_x[,1]+period*(counter-1)
}
sintro_1_x<-sinpks_1_x+period/2
sintro_2_x<-sinpks_2_x+period/2

le1<-5
le2<-8

#set up the panel for the Moran effects of the first env var
par(fig=c((yaxwd)/totwd,
          (yaxwd+panwd)/totwd,
          (xaxht+panht+gap+panht+titlespace+xaxht+panht+gap)/totht,
          (xaxht+panht+gap+panht+titlespace+xaxht+panht+gap+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)

#do the start of the plotting
ylimits<-range(envvar_1)
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
plot(tim,envvar_1[1,],type="n",ylim=ylimits,xaxs="i",xaxt="n",yaxt="n")
axis(side=1,label=FALSE)
#put the population rectangles, only done here so other stuff plots over it
xlimits<-range(tim)
rect(xlimits[1],vertoffset/2-1,xlimits[2],vertoffset/2+1,col="lightgrey",border=NA)
rect(xlimits[1],3*vertoffset/2-1,xlimits[2],3*vertoffset/2+1,col="lightgrey",border=NA)
rect(xlimits[1],5*vertoffset/2-1,xlimits[2],5*vertoffset/2+1,col="lightgrey",border=NA)
#put rectangles to highlight the synchrony, only done here so other stuff plots over it
ep<-1.2
rect(mean(sintro_1_x[,2])+le1-ep,-0.5,mean(sintro_1_x[,2])+le1+ep,2*vertoffset+4.5,border=NA,col=rgb(255/255,192/255,203/255,.5)) #"pink" is red    255, green  192, blue   203
rect(mean(sinpks_1_x[,3])+le1-ep,-0.5,mean(sinpks_1_x[,3])+le1+ep,2*vertoffset+4.5,border=NA,col=rgb(135/255,206/255,250/255,.5)) #"lightskyblue" is red    135, green  206, blue   250
#now finish plotting the env vars and proceed
lines(tim,envvar_1[1,],type="l")
lines(tim,envvar_1[2,],type="l")
lines(tim,envvar_1[3,],type="l")
text(10,vertoffset/2,latex2exp::TeX("population, $w_3(t)"),adj=c(0,.5),cex=0.6)
text(10,3*vertoffset/2,latex2exp::TeX("population, $w_2(t)"),adj=c(0,.5),cex=0.6)
text(10,5*vertoffset/2,latex2exp::TeX("population, $w_1(t)"),adj=c(0,.5),cex=0.6)

#add arrows showing max positive and negative effects on the pops
for (csw in 1:dim(sinpks_1_x)[1])
{
  for (cpk in 1:dim(sinpks_1_x)[2])
  {
    #plot the max positive effects on the pops as black arrows
    sx<-sinpks_1_x[csw,cpk]
    ex<-sx
    sy<-vertoffset*(csw-1)+1
    ey<-sy+.5
    lines(c(sx,ex),c(sy,ey))
    
    sx<-ex
    sy<-ey
    ex<-sx+le1
    ey<-sy
    lines(c(sx,ex),c(sy,ey))
    
    sx<-ex
    xy<-ey
    ex<-sx
    ey<-sy+1
    #lines(c(sx,ex),c(sy,ey))
    shape::Arrows(sx,sy,ex,ey,arr.length=0.15)
    
    #plot the max negative effects on the pops as red arrows
    sx<-sintro_1_x[csw,cpk]
    ex<-sx
    sy<-vertoffset*(csw-1)-1
    ey<-sy+1.75
    lines(c(sx,ex),c(sy,ey),col="red")
    
    sx<-ex
    sy<-ey
    ex<-sx+le1
    ey<-sy
    lines(c(sx,ex),c(sy,ey),col="red")
    
    sx<-ex
    xy<-ey
    ex<-sx
    ey<-sy+1.75
    #lines(c(sx,ex),c(sy,ey))
    shape::Arrows(sx,sy,ex,ey,arr.length=0.15,col="red")
  }
}

#label the arrows to show the lags
text(sinpks_1_x[3,2]+le1/2,2*vertoffset+1+.5,latex2exp::TeX("$l_{e1}$"),cex=0.7,adj=c(0.5,-.1))
text(sintro_1_x[2,3]+le1/2,vertoffset-1+1.75,latex2exp::TeX("$l_{e1}$"),cex=0.7,adj=c(0.5,-.1))
text(sinpks_1_x[1,4]+le1/2,1+.5,latex2exp::TeX("$l_{e1}$"),cex=0.7,adj=c(0.5,-.1))

#panel label and title
text(xlimits[1],ylimits[2],"(a)",adj=c(1.1,0),xpd=NA)
mtext("Synergistic Moran interactions",side=3,line=0.1,cex=1,at=xlimits[1],adj=c(-.1))

#label ln
par(xpd=NA)
lines(rep(sinpks_1_x[1,1],2),c(1,-4.5),type="l",lty="solid",col="violet")
par(xpd=FALSE)
mtext(latex2exp::TeX("$l_{n}$"),side=1,line=0,at=mean(c(sinpks_1_x[1,1],sinpks_2_x[3,1])),cex=0.7,adj=c(.5,.5))

#box around the whole thing
par(fig=c(0.005,0.995,.505,0.995),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(1:3,1:3,type="n",xaxt="n",yaxt="n")


#set up the panel for the Moran effects of the second env var
par(fig=c((yaxwd)/totwd,
          (yaxwd+panwd)/totwd,
          (xaxht+panht+gap+panht+titlespace+xaxht)/totht,
          (xaxht+panht+gap+panht+titlespace+xaxht+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

#do the start of the plotting
ylimits<-range(envvar_2)
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
plot(tim,envvar_2[1,],type="n",ylim=ylimits,lty="dashed",xaxs="i",xaxt="n",yaxt="n")
axis(side=1,label=TRUE)
mtext("Time",side=1,line=1.1)
#put the population rectangles, only done here so other stuff plots over it
xlimits<-range(tim)
rect(xlimits[1],vertoffset/2-1,xlimits[2],vertoffset/2+1,col="lightgrey",border=NA)
rect(xlimits[1],3*vertoffset/2-1,xlimits[2],3*vertoffset/2+1,col="lightgrey",border=NA)
rect(xlimits[1],5*vertoffset/2-1,xlimits[2],5*vertoffset/2+1,col="lightgrey",border=NA)
#put rectangles to highlight the synchrony, only done here so other stuff plots over it
ep<-1.2
rect(mean(sintro_2_x[,2])+le2-ep,-0.5,mean(sintro_2_x[,2])+le2+ep,2*vertoffset+4.5,border=NA,col=rgb(255/255,192/255,203/255,.5)) #"pink" is red    255, green  192, blue   203
rect(mean(sinpks_2_x[,3])+le2-ep,-0.5,mean(sinpks_2_x[,3])+le2+ep,2*vertoffset+4.5,border=NA,col=rgb(135/255,206/255,250/255,.5)) #"lightskyblue" is red    135, green  206, blue   250
#now finish plotting the env vars and proceed
lines(tim,envvar_2[1,],type="l",lty="dashed")
lines(tim,envvar_2[2,],type="l",lty="dashed")
lines(tim,envvar_2[3,],type="l",lty="dashed")
text(10,vertoffset/2,latex2exp::TeX("population, $w_3(t)"),adj=c(0,.5),cex=0.6)
text(10,3*vertoffset/2,latex2exp::TeX("population, $w_2(t)"),adj=c(0,.5),cex=0.6)
text(10,5*vertoffset/2,latex2exp::TeX("population, $w_1(t)"),adj=c(0,.5),cex=0.6)

for (csw in 1:dim(sinpks_2_x)[1])
{
  for (cpk in 1:dim(sinpks_2_x)[2])
  {
    #plot the max positive effects of the env var on the pops
    sx<-sinpks_2_x[csw,cpk]
    ex<-sinpks_2_x[csw,cpk]
    sy<-vertoffset*(csw-1)+1
    ey<-sy+.5
    lines(c(sx,ex),c(sy,ey))
    
    sx<-ex
    sy<-ey
    ex<-sx+le2
    ey<-sy
    lines(c(sx,ex),c(sy,ey))
    
    sx<-ex
    xy<-ey
    ex<-sx
    ey<-sy+1
    shape::Arrows(sx,sy,ex,ey,arr.length=0.15)
    
    #plot the max negative effects on the pops as red arrows
    sx<-sintro_2_x[csw,cpk]
    ex<-sx
    sy<-vertoffset*(csw-1)-1
    ey<-sy+2.2
    lines(c(sx,ex),c(sy,ey),col="red")
    
    sx<-ex
    sy<-ey
    ex<-sx+le2
    ey<-sy
    lines(c(sx,ex),c(sy,ey),col="red")
    
    sx<-ex
    xy<-ey
    ex<-sx
    ey<-sy+1.3
    #lines(c(sx,ex),c(sy,ey))
    shape::Arrows(sx,sy,ex,ey,arr.length=0.15,col="red")
  }
}

#label the arrows to show the lags
text(sinpks_2_x[3,2]+le2/2,2*vertoffset+1+.5,latex2exp::TeX("$l_{e2}$"),cex=0.7,adj=c(0.5,-.1))
text(sintro_2_x[2,3]+le2/2,vertoffset-1+2.2,latex2exp::TeX("$l_{e2}$"),cex=0.7,adj=c(0.5,-.1))
text(sinpks_2_x[1,4]+le2/2,1+.5,latex2exp::TeX("$l_{e2}$"),cex=0.7,adj=c(0.5,-.1))

#panel label
text(xlimits[1],ylimits[2],"(b)",adj=c(1.1,0),xpd=NA)

#label ln
par(xpd=NA)
lines(rep(sinpks_2_x[3,1],2),c(2*vertoffset+1+.5,2*vertoffset+1+7.5),type="l",lty="solid",col="violet")
par(xpd=FALSE)

#***part 2, which shows antagonistic effects, compressed presentation

#get all the quantities to plot
tim<-seq(from=0,to=80,length.out=1001)
period<-20

phase1<-(-1.5)
phase2<-phase1+2*pi*ln/period

phase_jitters_1<-2*pi*(2*runif(3)-1)/period
phase_jitters_2<-2*pi*(3*runif(3)-1)/period

vertoffset<-8

envvar_1<-matrix(NA,3,length(tim))
envvar_2<-matrix(NA,3,length(tim))
for (counter in 1:(dim(envvar_1)[1]))
{
  envvar_1[counter,]<-sin(2*pi*tim/period+phase1+phase_jitters_1[counter])+(counter-1)*vertoffset
  envvar_2[counter,]<-sin(2*pi*tim/period+phase2+phase_jitters_2[counter])+(counter-1)*vertoffset
}

sinpks_1_x<-matrix(NA,3,4)
sinpks_1_x[,1]<-(pi/2-phase1-phase_jitters_1)*period/2/pi
sinpks_2_x<-matrix(NA,3,4)
sinpks_2_x[,1]<-(pi/2-phase2-phase_jitters_2)*period/2/pi
for (counter in 2:4)
{
  sinpks_1_x[,counter]<-sinpks_1_x[,1]+period*(counter-1)
  sinpks_2_x[,counter]<-sinpks_2_x[,1]+period*(counter-1)
}
sintro_1_x<-sinpks_1_x+period/2
sintro_2_x<-sinpks_2_x+period/2

le1<-3
le2<-16

#set up the panel
par(fig=c((yaxwd)/totwd,
          (yaxwd+panwd)/totwd,
          (xaxht+panht+gap)/totht,
          (xaxht+panht+gap+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

#do the start of the plotting
ylimits<-range(envvar_1)
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
plot(tim,envvar_1[1,],type="n",ylim=ylimits,xaxs="i",xaxt="n",yaxt="n")
axis(side=1,label=FALSE)
#put the population rectangles, only done here so other stuff plots over it
xlimits<-range(tim)
rect(xlimits[1],vertoffset/2-1,xlimits[2],vertoffset/2+1,col="lightgrey",border=NA)
rect(xlimits[1],3*vertoffset/2-1,xlimits[2],3*vertoffset/2+1,col="lightgrey",border=NA)
rect(xlimits[1],5*vertoffset/2-1,xlimits[2],5*vertoffset/2+1,col="lightgrey",border=NA)
#put rectangles to highlight the synchrony, only done here so other stuff plots over it
ep<-1.2
rect(mean(sintro_1_x[,2])+le1-ep,-0.5,mean(sintro_1_x[,2])+le1+ep,2*vertoffset+4.5,border=NA,col=rgb(255/255,192/255,203/255,.5)) #"pink" is red    255, green  192, blue   203
rect(mean(sinpks_1_x[,3])+le1-ep,-0.5,mean(sinpks_1_x[,3])+le1+ep,2*vertoffset+4.5,border=NA,col=rgb(135/255,206/255,250/255,.5)) #"lightskyblue" is red    135, green  206, blue   250
#now finish plotting the env vars and proceed
lines(tim,envvar_1[1,],type="l")
lines(tim,envvar_1[2,],type="l")
lines(tim,envvar_1[3,],type="l")
text(10,vertoffset/2,latex2exp::TeX("population, $w_3(t)"),adj=c(0,.5),cex=0.6)
text(10,3*vertoffset/2,latex2exp::TeX("population, $w_2(t)"),adj=c(0,.5),cex=0.6)
text(10,5*vertoffset/2,latex2exp::TeX("population, $w_1(t)"),adj=c(0,.5),cex=0.6)

#add arrows showing max positive and negative effects on the pops
for (csw in 1:dim(sinpks_1_x)[1])
{
  for (cpk in 1:dim(sinpks_1_x)[2])
  {
    #plot the max positive effects on the pops as black arrows
    sx<-sinpks_1_x[csw,cpk]
    ex<-sx
    sy<-vertoffset*(csw-1)+1
    ey<-sy+.5
    lines(c(sx,ex),c(sy,ey))
    
    sx<-ex
    sy<-ey
    ex<-sx+le1
    ey<-sy
    lines(c(sx,ex),c(sy,ey))
    
    sx<-ex
    xy<-ey
    ex<-sx
    ey<-sy+1
    #lines(c(sx,ex),c(sy,ey))
    shape::Arrows(sx,sy,ex,ey,arr.length=0.15)
    
    #plot the max negative effects on the pops as red arrows
    sx<-sintro_1_x[csw,cpk]
    ex<-sx
    sy<-vertoffset*(csw-1)-1
    ey<-sy+1.75
    lines(c(sx,ex),c(sy,ey),col="red")
    
    sx<-ex
    sy<-ey
    ex<-sx+le1
    ey<-sy
    lines(c(sx,ex),c(sy,ey),col="red")
    
    sx<-ex
    xy<-ey
    ex<-sx
    ey<-sy+1.75
    #lines(c(sx,ex),c(sy,ey))
    shape::Arrows(sx,sy,ex,ey,arr.length=0.15,col="red")
  }
}

#label the arrows to show the lags
text(sinpks_1_x[3,2]+le1/4,2*vertoffset+1+.5,latex2exp::TeX("$l_{e1}$"),cex=0.7,adj=c(0.5,-.1))
text(sintro_1_x[2,3]+le1/4,vertoffset-1+1.75,latex2exp::TeX("$l_{e1}$"),cex=0.7,adj=c(0.5,-.1))
text(sinpks_1_x[1,4]+le1/4,1+.5,latex2exp::TeX("$l_{e1}$"),cex=0.7,adj=c(0.5,-.1))

#panel label and title
text(xlimits[1],ylimits[2],"(c)",adj=c(1.1,0),xpd=NA)
mtext("Antagonistic Moran interactions",side=3,line=0.1,cex=1,at=xlimits[1],adj=c(-.1))

#label ln
par(xpd=NA)
lines(rep(sinpks_1_x[1,1],2),c(1,-4.5),type="l",lty="solid",col="violet")
par(xpd=FALSE)
mtext(latex2exp::TeX("$l_{n}$"),side=1,line=0,at=mean(c(sinpks_1_x[1,1],sinpks_2_x[3,1])),cex=0.7,adj=c(.5,.5))

#box around the whole thing
par(fig=c(0.005,0.995,.005,0.495),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(1:3,1:3,type="n",xaxt="n",yaxt="n")

#set up the panel
par(fig=c((yaxwd)/totwd,
          (yaxwd+panwd)/totwd,
          (xaxht)/totht,
          (xaxht+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

#do the start of the plotting
ylimits<-range(envvar_2)
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
plot(tim,envvar_2[1,],type="n",ylim=ylimits,lty="dashed",xaxs="i",xaxt="n",yaxt="n")
axis(side=1,label=TRUE)
mtext("Time",side=1,line=1.1)
#put the population rectangles, only done here so other stuff plots over it
xlimits<-range(tim)
rect(xlimits[1],vertoffset/2-1,xlimits[2],vertoffset/2+1,col="lightgrey",border=NA)
rect(xlimits[1],3*vertoffset/2-1,xlimits[2],3*vertoffset/2+1,col="lightgrey",border=NA)
rect(xlimits[1],5*vertoffset/2-1,xlimits[2],5*vertoffset/2+1,col="lightgrey",border=NA)
#put rectangles to highlight the synchrony, only done here so other stuff plots over it
ep<-1.2
rect(mean(sintro_2_x[,2])+le2-ep,-0.5,mean(sintro_2_x[,2])+le2+ep,2*vertoffset+4.5,border=NA,col=rgb(255/255,192/255,203/255,.5)) #"pink" is red    255, green  192, blue   203
rect(mean(sinpks_2_x[,2])+le2-ep,-0.5,mean(sinpks_2_x[,2])+le2+ep,2*vertoffset+4.5,border=NA,col=rgb(135/255,206/255,250/255,.5)) #"lightskyblue" is red    135, green  206, blue   250
#now finish plotting the env vars and proceed
lines(tim,envvar_2[1,],type="l",lty="dashed")
lines(tim,envvar_2[2,],type="l",lty="dashed")
lines(tim,envvar_2[3,],type="l",lty="dashed")
text(10,vertoffset/2,latex2exp::TeX("population, $w_3(t)"),adj=c(0,.5),cex=0.6)
text(10,3*vertoffset/2,latex2exp::TeX("population, $w_2(t)"),adj=c(0,.5),cex=0.6)
text(10,5*vertoffset/2,latex2exp::TeX("population, $w_1(t)"),adj=c(0,.5),cex=0.6)

for (csw in 1:dim(sinpks_2_x)[1])
{
  for (cpk in 1:dim(sinpks_2_x)[2])
  {
    #plot the max positive effects of the env var on the pops
    sx<-sinpks_2_x[csw,cpk]
    ex<-sinpks_2_x[csw,cpk]
    sy<-vertoffset*(csw-1)+1
    ey<-sy+.5
    lines(c(sx,ex),c(sy,ey))
    
    sx<-ex
    sy<-ey
    ex<-sx+le2
    ey<-sy
    lines(c(sx,ex),c(sy,ey))
    
    sx<-ex
    xy<-ey
    ex<-sx
    ey<-sy+1
    shape::Arrows(sx,sy,ex,ey,arr.length=0.15)
    
    #plot the max negative effects on the pops as red arrows
    sx<-sintro_2_x[csw,cpk]
    ex<-sx
    sy<-vertoffset*(csw-1)-1
    ey<-sy+2.2
    lines(c(sx,ex),c(sy,ey),col="red")
    
    sx<-ex
    sy<-ey
    ex<-sx+le2
    ey<-sy
    lines(c(sx,ex),c(sy,ey),col="red")
    
    sx<-ex
    xy<-ey
    ex<-sx
    ey<-sy+1.3
    #lines(c(sx,ex),c(sy,ey))
    shape::Arrows(sx,sy,ex,ey,arr.length=0.15,col="red")
  }
}

#label the arrows to show the lags
text(sinpks_2_x[3,2]+3*le2/5,2*vertoffset+1+.5,latex2exp::TeX("$l_{e2}$"),cex=0.7,adj=c(0.5,-.1))
#text(sintro_2_x[2,3]+le2/2,vertoffset-1+2.2,latex2exp::TeX("$l_{e2}$"),cex=0.7,adj=c(0.5,-.1))
#text(sinpks_2_x[1,4]+le2/2,1+.5,latex2exp::TeX("$l_{e2}$"),cex=0.7,adj=c(0.5,-.1))

#panel label
text(xlimits[1],ylimits[2],"(d)",adj=c(1.1,0),xpd=NA)

#label ln
par(xpd=NA)
lines(rep(sinpks_2_x[3,1],2),c(2*vertoffset+1+.5,2*vertoffset+1+7.5),type="l",lty="solid",col="violet")
par(xpd=FALSE)

dev.off()
