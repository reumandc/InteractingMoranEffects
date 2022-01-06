#This script implements very simple theoretical model explored in the Introduction

rm(list=ls())

#***
#External codes needed
#***

#packages needed (invoked with "::"): matrixcalc, shape, latex2exp

source("IntroTheoryIA.R")

#***
#Location for storing results and other prep
#***

resloc<-"../Results/Theory/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***
#Use the function to explore some examples of interactions, looking for a synergistic example
#and a destuctive example, and save those results for possible inclusion in the Intro
#***

a<-1
b<-1

#synergistic example
covmat<-matrix(c(1,.7,.5,.5,
                 .7,1,.5,.5,
                 .5,.5,1,.7,
                 .5,.5,.7,1),4,4,byrow=TRUE)
res_syn<-list(a=a,b=b,covmat=covmat,res=IntroTheoryIA(a,b,covmat))
saveRDS(res_syn,file=paste0(resloc,"IntroTheory_SynergisticExample.Rds"))

#destructive example
covmat<-matrix(c(1,.7,-.5,-.5,
                 .7,1,-.5,-.5,
                 -.5,-.5,1,.7,
                 -.5,-.5,.7,1),4,4,byrow=TRUE)
res_des<-list(a=a,b=b,covmat=covmat,res=IntroTheoryIA(a,b,covmat))
saveRDS(res_des,file=paste0(resloc,"IntroTheory_DestructiveExample.Rds"))

#***
#Now explore a one parameter family of possibilities and plot it
#***

#get the result
lo<-101
od<-seq(from=-1,to=1,length.out=lo)
res<-data.frame(od=od,sync=NA,sync.in=NA)
for (counter in 1:(dim(res)[1]))
{
  covmat<-matrix(c(1,.7,od[counter],od[counter],
                   .7,1,od[counter],od[counter],
                   od[counter],od[counter],1,.7,
                   od[counter],od[counter],.7,1),4,4,byrow=TRUE)
  h<-IntroTheoryIA(a=1,b=1,covmat=covmat)
  res[counter,2]<-h['sync']
  res[counter,3]<-h['sync.ni']
}

#dimensions for the figure to be generated, inches
totwd<-3.5
totht<-totwd
xaxht<-.5
yaxwd<-.5
gap<-0.1
panwd<-(totwd-gap-yaxwd)
panht<-(totht-gap-xaxht)

pdf(file=paste0(resloc,"IntroTheory.pdf"),height=totht,width=totwd)

par(fig=c((yaxwd)/totwd,
          (yaxwd+panwd)/totwd,
          (xaxht)/totht,
          (xaxht+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)

plot(res$od,res$sync,type="l",xlab="Parameter",ylab="Population correlation")
lines(res$od,res$sync.in,type="l",lty="dashed")

mtext(latex2exp::TeX("Parameter, $\\phi$"),1,1.2)
mtext("Population synchrony (correlation)",2,1.2)

inds<-which(res$sync<=res$sync.in)
x<-c(res$od[inds],rev(res$od[inds]))
y<-c(res$sync.in[inds],rev(res$sync[inds]))
polygon(x,y,col="pink",border=NA)

inds<-which(res$sync>=res$sync.in)
x<-c(res$od[inds],rev(res$od[inds]))
y<-c(res$sync[inds],rev(res$sync.in[inds]))
polygon(x,y,col="lightskyblue",border=NA)

lines(res$od,res$sync,type="l")
lines(res$od,res$sync.in,type="l",lty="dashed")

shape::Arrows(-.25,0,-.6,.5)
text(-.25,0,"Negative \n interaction \n effects",adj=c(0.5,1))
shape::Arrows(.5,.5,.7,.72)
text(.5,.5,"Positive \n interaction \n effects",adj=c(0.5,1))

text(-1,.75,"(a)",cex=1.1,adj=c(0,0))

legend(x="bottomright",legend=c("Pop. sync.","Pop. sync. with interaction \n effects removed"),cex=.85,lty=c("solid","dashed"))

dev.off()