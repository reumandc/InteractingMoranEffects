#A function for computing the mean Pearson correlation between rows of two matrices, with a lag
#
#Args
#m1, m2       The two matrices, must be the same size
#lag          The lag. Positive values only. Means you're computing correlations between m1[i,(lag+1):(dim(m1)[2])] and m2[i,1:(dim(m2)[2]-lag)]
#
#Output
#A single number, the mean correlation between rows of m1 and corresponding rows of m2
mncor_tiny<-function(m1,m2,lag)
{
  numts<-dim(m1)[1]
  lents<-dim(m1)[2]

  #error catching
  if (any(dim(m1)!=dim(m2)))
  {
    stop("Error in mncor_tiny: m1 and m2 must have the same dimension")
  }
  if (lag<0)
  {
    stop("Error in mncor_tiny: lag must be >=0")
  }
  
  res<-NA*numeric(numts)
  for (counter in 1:numts)
  {
    res[counter]<-cor(m1[counter,(lag+1):lents],m2[counter,1:(lents-lag)],method="pearson")
  }
  
  return(mean(res))
}