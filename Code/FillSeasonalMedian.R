#A function for filling NAs in a quarterly time series with the seasonal median
#
#Args
#x      A quarterly time series
#
#Output - the same time series, but with NAs filled
#
fill_with_seasonal_median<-function(x)
{
  allinds<-1:length(x)
  nainds<-which(is.na(x))
  res<-x
  for (counter in 1:length(nainds))
  {
    res[nainds[counter]]<-median(x[allinds[(allinds %% 4)==(nainds[counter] %% 4)]],na.rm=TRUE)
  }
  return(res)
}
