#This tests the function in PlotPhaseFunc.R, but since that's a plotting function this script is
#not covered by the automated testing through the testthat package, which is why the file name
#does not start with "tests".

source("PlotPhaseFunc.R")

#a simple case not ever going through phase pi=-pi
t<-1:10
z<-complex(modulus=1:10,argument=seq(from=-3,to=3,length.out=10))
plot(0,0,xlim=range(t),ylim=c(-pi,pi))
plotphasefunc(t,z)

#a case where you pass through phase pi=-pi from the pi side to the -pi side
z<-complex(modulus=1:5,argument=c(pi*seq(from=.1,to=.9,length.out=10),-pi*seq(from=.9,to=.1,length.out=10)))
t<-1:length(z)
plot(0,0,xlim=range(t),ylim=c(-pi,pi))
plotphasefunc(t,z)

# case where you pass through phase pi=-pi from the -pi side to the pi side
z<-complex(modulus=1:5,argument=c(-pi*seq(from=.1,to=.9,length.out=10),pi*seq(from=.9,to=.1,length.out=10)))
t<-1:length(z)
plot(0,0,xlim=range(t),ylim=c(-pi,pi))
plotphasefunc(t,z)

#a random case
set.seed(102)
t<-1:20
z<-complex(modulus=1:20,argument=runif(20,min=-pi,max=pi))
plot(0,0,xlim=range(t),ylim=c(-pi,pi))
plotphasefunc(t,z)

