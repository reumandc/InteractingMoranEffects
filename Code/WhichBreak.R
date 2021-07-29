#Takes a set of breaks as in "hist" and a set of values and returns a set of indices of 
#the same length as vals indicating between which two breaks each value is.
#
#Args
#breaks         A set of breaks, as for the function hist - evenly spaced increasing values
#vals           A set of values, assumed to be >breaks[1] and <=breaks[length(breaks)]
#
#Output
#A vector of the same length as vals with numbers between 1 and length(breaks)-1 indicating
#which of the breaks the corresponding value fell between.
#
WhichBreak<-function(breaks,vals)
{
  if (!isTRUE(all.equal(diff(breaks),rep(breaks[2]-breaks[1],length(breaks)-1))))
  {
    stop("Error in WhichBreak: breaks must have evenly spaced values")
  }
  if (breaks[2]-breaks[1]<=0)
  {
    stop("Error in WhichBreak: breaks must have increasing values")
  }
  if (any(vals<=breaks[1]) || any(vals>breaks[length(breaks)]))
  {
    stop("Error in WhichBreak: some entries of vals fall outside the range of breaks")
  }
  
  res<-NA*numeric(length(vals))
  for (counter in 1:length(vals))
  {
    res[counter]<-max(which(breaks<vals[counter]))
  }
  
  return(res)
}