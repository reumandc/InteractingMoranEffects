#Computes the mean coefficient of determination of a linear regression computed for each location, as well as
#mean coefficients.
#
#Args
#kelp           Locations by time matrix for kelp data
#NO3            Same for nitrates
#waves          Same for waves
#
#Output
#For each location, does the linear regression with lags 0 and 1 for waves and NO3, and gets the coefficient of 
#determination (R^2). Returns the average of these values across all locations. Also returns the average NO3 
#coefficients in the regressions, also returns the average waves coefficients in the regressions.
#
mncoefdet_tiny2<-function(kelp,NO3,waves)
{
  if (any(dim(kelp)!=dim(NO3)) || any(dim(kelp)!=dim(waves)))
  {
    stop("Error in mncoefdet_tiny2: kelp, NO3 and waves input matrices must be the same dimensions")
  }

  numts<-dim(kelp)[1]
  lents<-dim(kelp)[2]
  
  res_Rsq<-NA*numeric(numts)
  res_NO3coefs<-matrix(NA,numts,2)
  res_wavescoefs<-matrix(NA,numts,2)
  for (counter in 1:numts)
  {
    #extract time series for this location
    tskelp<-kelp[counter,]
    tsNO3<-NO3[counter,]
    tswaves<-waves[counter,]
    
    #implement the desired lags
    tskelp<-tskelp[2:lents]
    tsNO3_l0<-tsNO3[2:lents]
    tsNO3_l1<-tsNO3[1:(lents-1)]
    tswaves_l0<-tswaves[2:lents]
    tswaves_l1<-tswaves[1:(lents-1)]
    
    #do the regressiom
    mod<-lm(tskelp~tsNO3_l0+tsNO3_l1+tswaves_l0+tswaves_l1)
    res_Rsq[counter]<-summary(mod)$r.squared
    res_NO3coefs[counter,]<-unname(coef(mod)[2:3])
    res_wavescoefs[counter,]<-unname(coef(mod)[4:5])
  }
  
  return(c(mean(res_Rsq),apply(FUN=mean,X=res_NO3coefs,MARGIN=2),apply(FUN=mean,X=res_wavescoefs,MARGIN=2)))
}
