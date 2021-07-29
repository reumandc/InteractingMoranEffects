#This tests the function in plotColLine.R, but since that's a plotting function this script is
#not covered by the automated testing through the testthat package, which is why the file name
#does not start with "tests".

source("plotColLine.R")

x<-1:5
y<-c(0,1,0,1,0)
cols<-c("red","blue","green","yellow","purple")
plotColLine(x,y,cols,lnwid=2)
