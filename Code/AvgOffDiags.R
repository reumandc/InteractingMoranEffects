#Takes an array of dimensions n by n by m and averages the off-diagonal vectors.
#
#Args
#A          The array
#
#Output
#A vector of length m
#
avgoffdiags<-function(A)
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
  
  n<-dim(A)[1]
  m<-dim(A)[3]
  for (counter in 1:m)
  {
    h<-A[,,counter]
    diag(h)<-NA
    A[,,counter]<-h
  }
  
  res<-apply(FUN=mean,X=A,MARGIN=3,na.rm=TRUE)
  return(res)
}