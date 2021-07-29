#Plots a line, but with shifting colors.
#
#Args
#x, y       Vectors of the same length to plot
#cols       Same length as x and y, colors to use
#lnwid      Line width to use      
#...        Passed to "plot" command, used for controlling axis limits and the like
#
#Notes
#Makes a line plot using the colors on the corresponding parts of the plot.
#
plotColLine<-function(x,y,cols,lnwid,...)
{
  lx<-length(x)
  
  plot(x,y,type="n",...)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey90")
  
  lines(c(x[1],mean(x[1:2])),c(y[1],mean(y[1:2])),lwd=lnwid,col=cols[1])
  for (counter in 2:(length(x)-1))
  {
    lines(c(mean(c(x[counter-1],x[counter])),x[counter]),
          c(mean(c(y[counter-1],y[counter])),y[counter]),lwd=lnwid,col=cols[counter])
    lines(c(x[counter],mean(c(x[counter],x[counter+1]))),
          c(y[counter],mean(c(y[counter],y[counter+1]))),lwd=lnwid,col=cols[counter])
  }
  lines(c(mean(x[(lx-1):lx]),x[lx]),c(mean(y[(lx-1):lx]),y[lx]),lwd=lnwid,col=cols[lx])
}

