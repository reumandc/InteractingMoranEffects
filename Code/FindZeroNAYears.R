#A function for finding how many years of data are comprised solely of NA/0 values.
#
#Args
#x          A time series, assumed quarterly, assumed to start with q1 of a year and end with q4 of a year
#
#Output
#A single number which is the number of calendar years for which data were comprised solely of 0/NA values
#
find_zero_NA_years<-function(x)
{
  if (length(x)%%4 != 0)
  {
    stop("Error in find_zero_NA_years: x needs to have length divisible by four")
  }
  
  numyrs<-0
  for (counter in 1:(length(x)/4))
  {
    h<-x[(4*counter-3):(4*counter)]
    h<-h[!is.na(h)]
    if (all(h==0))
    {
      numyrs<-numyrs+1
    }
  }
  
  return(numyrs)
}