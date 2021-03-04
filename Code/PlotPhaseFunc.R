#This function plots the phase of a complex-valued function of a real parameter on an axis already 
#assumed to exist.
#
#Args
#t          The real parameter. Assumed to be strictly increasing values.
#z          Complex numbers corresponding to the entries of t
#pch, cex   Plotting parameters for the points
#lty        Line type for the connecting lines
#col        Color for plotting
#
#Output
#The line connecting between Arg(z[a]) and Arg(z[a+1]) is assumed to go the shorter of the two routes
#around the phase circle. In the event z[a] and z[a+1] are 180 degrees apart, the route through 1
#on the complex plane is used.
#
plotphasefunc<-function(t,z,pch=20,cex=0.5,lty="solid",col="black")
{
  #error catching
  if (class(z)!="complex")
  {
    stop("Error in plotphasefunc: z must be a complex vector")
  }
  if (class(t)!="numeric" && class(t)!="integer")
  {
    stop("Error in plotphasefunc: t must be a numeric vector")
  }
  if (length(t)!=length(z))
  {
    stop("Error in plotphasefunc: t and z arguments must be vectors of the same length")
  }
  if (any(diff(t)<=0))
  {
    stop("Error in plotphasefunc: values of t must be strictly increasing")
  }
  
  #plot the points (without the lines connecting them)
  points(t,Arg(z),type="p",pch=pch,cex=cex,col=col)
  
  #now plot the lines between them
  for (counter in 1:(length(t)-1))
  {
    a1<-Arg(z[counter])
    a2<-Arg(z[counter+1])
    t1<-t[counter]
    t2<-t[counter+1]
    if (abs(a1-a2)<=pi)
    {
      #this is the case where you pass through phase 0
      lines(c(t1,t2),c(a1,a2),type="l",lty=lty,col=col)      
    }else
    {
      #this is the case where you pass through phase pi=-pi
      if (a1>0 && a2<0)
      {
        #this is where you pass through pi=-pi from the pi side to the -pi side
        d_a1tom1<-pi-a1
        d_m1toa2<-a2-(-pi)
        tm<-(d_a1tom1/(d_a1tom1+d_m1toa2))*(t2-t1)+t1
        lines(c(t1,tm),c(a1,pi),type="l",lty=lty,col=col)
        lines(c(tm,t2),c(-pi,a2),type="l",lty=lty,col=col)
        next
      }
      if (a1<0 && a2>0)
      {
        #this is where you pass through pi=-pi from the -pi side to the pi side
        d_a1tom1<-a1-(-pi)
        d_m1toa2<-pi-a2
        tm<-(d_a1tom1/(d_a1tom1+d_m1toa2))*(t2-t1)+t1
        lines(c(t1,tm),c(a1,-pi),type="l",lty=lty,col=col)
        lines(c(tm,t2),c(pi,a2),type="l",lty=lty,col=col)
        next
      }
      stop("Error in plotphasefunc: the code reached a place I thought it would never reach")
    }
  }
  
  return(NULL)
}