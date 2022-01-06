#This function calculates the strength of interaction effects in the simple example of SI section S1
#of Sheppard et al (2019, Plos Comp Biol). Notation is the same as that paper.
#
#Args
#a, b       The coefficients, see Sheppard et al (2019)
#covmat     The covariance matrix of (alpha(1),alpha(2),beta(1),beta(2)). Must have 1s down the 
#             diagonal, if only for consistency with the same simplifying assumption made in
#             Sheppard et al (2019), SI section S1.
#
#Output - a length-2 named vector with these components 
#sync       The population synchrony of the model (measured with correlation)
#sync.ni    What the population synchrony would be if the two Moran drivers were independent.
#
#Note
#Returns NA for sync if covmat is not positive definite, returns NA for sync.in if the analogous
#covariance matrix for is is not positive definite.
#
IntroTheoryIA<-function(a,b,covmat)
{
  if (any(diag(covmat)!=1))
  {
    stop("Error in IntroTheoryIA: covmat should have diagonal entries equal to 1")
  }

  if (matrixcalc::is.positive.definite(covmat))
  {
    res1<-(a^2*covmat[1,2]+a*b*covmat[1,4]+a*b*covmat[2,3]+b^2*covmat[3,4])/
      (sqrt(a^2+b^2+2*a*b*covmat[1,3])*sqrt(a^2+b^2+2*a*b*covmat[2,4]))
  } else
  {
    res1<-NA  
  }
  
  redcovmat<-covmat
  redcovmat[1,3:4]<-0
  redcovmat[2,3:4]<-0
  redcovmat[3,1:2]<-0
  redcovmat[4,1:2]<-0
  if (matrixcalc::is.positive.definite(redcovmat))
  {
    res2<-(a^2*redcovmat[1,2]+a*b*redcovmat[1,4]+a*b*redcovmat[2,3]+b^2*redcovmat[3,4])/
      (sqrt(a^2+b^2+2*a*b*redcovmat[1,3])*sqrt(a^2+b^2+2*a*b*redcovmat[2,4]))
  } else
  {
    res2<-NA
  }
  
  return(c(sync=res1,sync.ni=res2))
}