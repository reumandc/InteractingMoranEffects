#Takes an array of dimensions n by n by m and averages the off-diagonal vectors. Prior to averaging,
#if the second argument is not NA, elements of the first argument are rescaled in a way that uses
#the second argument.
#
#Args
#A          The array
#resc       Vector of length n, for a certain kind of rescaling prior to averaging. Default
#             NA means don't perform any rescaling.
#
#Output
#A vector of length m
#
#Details
#Picture that A is Sww, or some component of that, and we are rescaling by geometric
#means of variances of time series. That's what this is for. The argument resc should
#just be the variances of the time series.
#
avgoffdiags<-function(A,resc=NA)
{
  #error catching
  if (dim(A)[1]!=dim(A)[2])
  {
    stop("Error in avgoffdiags: A must have dimensions n by n by m")
  }
  if (any(!is.finite(A)))
  {
    stop("Error in avgoffdiags: A must have all finite entries")
  }
  if (!is.na(resc[1]) && (length(resc)!=dim(A)[1] || !all(is.finite(resc))))
  {
    stop("Error in avgoffdiags: bad value for resc")
  }
  
  #get ready for the normalization, if desired
  if (!is.na(resc[1]))
  {
    rescsqrt<-sqrt(resc)
    rescmat<-outer(X=rescsqrt,Y=rescsqrt,FUN="*")
  }else
  {
    rescmat<-matrix(1,nrow(A),ncol(A))  
  }
  
  n<-dim(A)[1]
  m<-dim(A)[3]
  for (counter in 1:m)
  {
    h<-A[,,counter]
    diag(h)<-NA
    A[,,counter]<-h/rescmat
  }
  
  res<-apply(FUN=mean,X=A,MARGIN=3,na.rm=TRUE)
  return(res)
}