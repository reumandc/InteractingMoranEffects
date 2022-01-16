#This script makes a pedagogical figure. No modelling here, just presentation.

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): shape

#***
#Location for storing results and other prep
#***

resloc<-"../Results/Theory/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***
#Make the figure
#***

#dimensions for the figure to be generated, inches
totwd<-3.5
totht<-totwd
gap<-0.15
xaxht<-0.5
panwd<-(totwd-2*gap)
panht<-(totht-2*gap-xaxht)/2

pdf(file=paste0(resloc,"PedagogFig.pdf"),height=totht,width=totwd)

#***part 1, which shows synergistic effects

#set up the panel
par(fig=c((gap)/totwd,
          (gap+panwd)/totwd,
          (xaxht+gap+panht)/totht,
          (xaxht+gap+2*panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)

#get and plot the sin waves
tim<-seq(from=0,to=100,length.out=10001)
period<-20
phi1<-(-1.5)
phi2<-phi1-pi/2
sin1<-sin(2*pi*tim/period+phi1)
sin2<-sin(2*pi*tim/period+phi2)
offset1<-0
offset2<-5
y1<-sin1+offset1
y2<-sin2+offset2

xlimits<-range(tim)
ylimits<-range(y1,y2)
dy<-diff(ylimits)
ylimits[2]<-ylimits[2]+.15*dy
ylimits[1]<-ylimits[1]-.01*dy

plot(tim,y1,type="l",xlim=xlimits,ylim=ylimits,xaxs="i",xaxt="n",yaxt="n")
lines(tim,y2,type="l")

#label populations
ep<-0.4
ctr<-mean(c(offset1,offset2))
rect(xlimits[1],ctr-ep,xlimits[2],ctr+ep,col="lightgrey",border=NA)
text(17,ctr,latex2exp::TeX("population, $w_i(t)"),adj=c(0,.5),cex=0.8)

#illustrate the lag between the two environmental variables
lines(rep((pi/2-phi2)*period/2/pi,2),c(6,3.2),type="l",lty="dashed")
lines(rep((pi/2-phi1)*period/2/pi,2),c(1,3.8),type="l",lty="dashed")
text(mean(c((pi/2-phi2)*period/2/pi,(pi/2-phi1)*period/2/pi)),3.5,latex2exp::TeX("$l_n$"),adj=c(0.5,0.5),cex=0.8)

#illustrate the delay of env effects on population
mutdel<-2
lines(rep((pi/2-phi2)*period/2/pi+3*period,2),c(6,3.5),type="l",lty="solid")
lines(c((pi/2-phi2)*period/2/pi+3*period,(pi/2-phi2)*period/2/pi+3*period+period/4+mutdel),rep(3.5,2))
shape::Arrows((pi/2-phi2)*period/2/pi+3*period+period/4+mutdel,3.5,
              (pi/2-phi2)*period/2/pi+3*period+period/4+mutdel,ctr+.2,
              arr.length=.2)
text(.55*((pi/2-phi2)*period/2/pi+3*period+.75)+.45*((pi/2-phi2)*period/2/pi+3*period+period/4+mutdel),
      3.5,latex2exp::TeX("$l_{e1}$"),cex=0.8,adj=c(0.5,-.2))

lines(rep((pi/2-phi1)*period/2/pi+3*period,2),c(1,1.25),type="l",lty="solid")
lines(c((pi/2-phi1)*period/2/pi+3*period,(pi/2-phi1)*period/2/pi+3*period+mutdel+period/2),rep(1.25,2))
shape::Arrows((pi/2-phi1)*period/2/pi+3*period+mutdel+period/2,1.25,
              (pi/2-phi1)*period/2/pi+3*period+mutdel+period/2,ctr-.2,
              arr.length=.2)
text(.5*((pi/2-phi1)*period/2/pi+3*period)+.5*((pi/2-phi1)*period/2/pi+3*period+mutdel+period/2),
     1.25,latex2exp::TeX("$l_{e2}"),cex=0.8,adj=c(0.5,-.2))

#illustrate which way time goes
#shape::Arrows(20,6.8,90,6.8)
#text(20,6.8,latex2exp::TeX("time, $t$"),adj=c(0,-.2),cex=0.8)

#panel label
text(0,6,"(b)",cex=1.1,adj=c(-.2,0))

#***part 2, which shows antagonistic effects

#set up the panel
par(fig=c((gap)/totwd,
          (gap+panwd)/totwd,
          (xaxht)/totht,
          (xaxht+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

#get and plot the sin waves
tim<-seq(from=0,to=100,length.out=10001)
period<-20
phi1<-(-1.5)
phi2<-phi1-pi/2
sin1<-sin(2*pi*tim/period+phi1)
sin2<-sin(2*pi*tim/period+phi2)
offset1<-0
offset2<-5
y1<-sin1+offset1
y2<-sin2+offset2

xlimits<-range(tim)
ylimits<-range(y1,y2)
dy<-diff(ylimits)
ylimits[2]<-ylimits[2]+.15*dy
ylimits[1]<-ylimits[1]-.01*dy

plot(tim,y1,type="l",xlim=xlimits,ylim=ylimits,xaxs="i",yaxt="n")
lines(tim,y2,type="l")
mtext("Time step",1,1.2)

#label populations
ep<-0.4
ctr<-mean(c(offset1,offset2))
rect(xlimits[1],ctr-ep,xlimits[2],ctr+ep,col="lightgrey",border=NA)
text(17,ctr,latex2exp::TeX("population, $w_i(t)"),adj=c(0,.5),cex=0.8)

#illustrate the lag between the two environmental variables
lines(rep((pi/2-phi2)*period/2/pi,2),c(6,3.2),type="l",lty="dashed")
lines(rep((pi/2-phi1)*period/2/pi,2),c(1,3.8),type="l",lty="dashed")
text(mean(c((pi/2-phi2)*period/2/pi,(pi/2-phi1)*period/2/pi)),3.5,latex2exp::TeX("$l_n$"),adj=c(0.5,0.5),cex=0.8)

#illustrate the delay of env effects on population
mutdel<-17
lines(rep((pi/2-phi2)*period/2/pi+3*period,2),c(6,3.5),type="l",lty="solid")
lines(c((pi/2-phi2)*period/2/pi+3*period,(pi/2-phi2)*period/2/pi+3*period-period/4+mutdel),rep(3.5,2))
shape::Arrows((pi/2-phi2)*period/2/pi+3*period-period/4+mutdel,3.5,
              (pi/2-phi2)*period/2/pi+3*period-period/4+mutdel,ctr+.2,
              arr.length=.2)
text(.75*((pi/2-phi2)*period/2/pi+3*period+.75)+.25*((pi/2-phi2)*period/2/pi+3*period+period/4+mutdel),
     3.5,latex2exp::TeX("$l_{e1}$"),cex=0.8,adj=c(0.5,-.2))

lines(rep((pi/2-phi1)*period/2/pi+3*period,2),c(1,1.25),type="l",lty="solid")
lines(c((pi/2-phi1)*period/2/pi+3*period,(pi/2-phi1)*period/2/pi+3*period+mutdel-period/2),rep(1.25,2))
shape::Arrows((pi/2-phi1)*period/2/pi+3*period+mutdel-period/2,1.25,
              (pi/2-phi1)*period/2/pi+3*period+mutdel-period/2,ctr-.2,
              arr.length=.2)
text(.5*((pi/2-phi1)*period/2/pi+3*period)+.5*((pi/2-phi1)*period/2/pi+3*period+mutdel-period/2),
     1.25,latex2exp::TeX("$l_{e2}"),cex=0.8,adj=c(0.5,-.2))

#illustrate which way time goes
#shape::Arrows(20,6.8,90,6.8)
#text(20,6.8,latex2exp::TeX("time, $t$"),adj=c(0,-.2),cex=0.8)

#panel label
text(0,6,"(c)",cex=1.1,adj=c(-.2,0))

dev.off()
