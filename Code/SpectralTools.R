#Some tools for Fourier analysis

#A raw, unsmoothed periodogram. Detrending is optionally done by this function, the zero frequency is optionally
#removed, and the symmetric part of the result is optionally cut. 
#
#Args
#x            A time series as a vector
#detrend      TRUE or FALSE according to whether you want to linearly detrend the time series first
#removezero   TRUE or FALSE according to whether you want to remove the zero frequency from the result or not
#cutsym       TRUE or FALSE depending on whether you want to cut out the symmetric part of the periodogram
#
#Note: Normalization is done so that the mean across frequencies (using all nonzero canonical frequencies from 0 
#to 1) of the spectrum is the variance of the original time series. That's under the options detrend=F, 
#removezero=T, cutsym=F.
#
myspecraw<-function(x,detrend=TRUE,removezero=TRUE,cutsym=TRUE)
{
  lx<-length(x)

  if (detrend==TRUE)
  {
    tforx<-1:lx
    x<-stats::residuals(stats::lm(x~tforx))
  }

  h<-stats::fft(x)
  h<-Re(h*Conj(h))
  h<-h/lx
  freq<-(0:(lx-1))/lx

  if (removezero==TRUE)
  {
    h<-h[-1]
    freq<-freq[-1]
  }
  
  if (cutsym==TRUE)
  {
    h<-h[freq<=0.5]
    freq<-freq[freq<=0.5]
  }
  
  return(list(freq=freq,spec=h))
}

#A raw, unsmoothed periodogram estimator of the spectral matrix. Detrending is optionally done by this function, 
#the zero frequency is optionally removed, and the symmetric part of the result is optionally cut. 
#
#Args
#x            A multivariate time series as a matrix, each row containing a component time series
#detrend      TRUE or FALSE according to whether you want to linearly detrend each component time series first
#removezero   TRUE or FALSE according to whether you want to remove the zero frequency from the result or not
#cutsym       TRUE or FALSE depending on whether you want to cut out the symmetric part of the result, thereby
#               returning only the results for frequency <= 0.5.
#
#Note: Normalization is done so that the mean across frequencies (using all nonzero canonical frequencies from 0 
#to 1) of a spectrum/cospectrum is the variance/covariance of the original time series. That's under the options 
#detrend=F, removezero=T, cutsym=F.
#
myspecmatraw<-function(x,detrend=TRUE,removezero=TRUE,cutsym=TRUE)
{
  Tx<-dim(x)[2] #length of time series
  N<-dim(x)[1] #number of components
  
  #detrend
  if (detrend==TRUE)
  {
    tforx<-1:Tx
    for (counter in 1:N)
    {
      y<-x[counter,]
      x[counter,]<-stats::residuals(stats::lm(y~tforx))  
    }
  }
  
  #get the raw periodogram with the zero frequency removed 
  fftx<-mymvfft(x)
  I<-array(complex(real=0,imaginary=0),dim=c(N,N,Tx))
  for (a in 1:N)
  {
    for (b in 1:N)
    {
      I[a,b,]<-fftx[a,]*Conj(fftx[b,])
    }
  }
  
  if (removezero==TRUE)
  {
    I<-I[,,2:Tx]/Tx
    freq<-(1:(Tx-1))/Tx
  }else
  {
    I<-I/Tx
    freq<-(0:(Tx-1))/Tx
  }
  
  #cut the symmetric part of the spectrum, if desired
  if (cutsym==TRUE)
  {
    I<-I[,,freq<=0.5]
    freq<-freq[freq<=0.5]
  }
  
  return(list(freq=freq,spec=I))
}

#Takes the fft of each row of the matrix x. Just a convenience function.
#
#Args
#x      A matrix
#
mymvfft<-function(x)
{
  return(t(stats::mvfft(t(x))))
}

#Inverse of the above
#
#Args
#x      A matrix
#
imymvfft<-function(x)
{
  return(t(stats::mvfft(t(x),inverse=TRUE))/dim(x)[2])
}

#This function does a kind of smoothing on x, but with some wrinkles which are specific to the application to the 
#Brillinger consistent estimator of the spectrum.
#
#Args
#sm     A numeric vector of odd length representing the smoothing kernel. The middle one is the center of the kernel. 
#         Must have length less than that of x. 
#x      A vector to be smoothed. This is intended to be a raw periodogram, with the symmetric part retained but the 
#         zero frequency removed (removed, not set to 0).
#
#Output - the smoothed periodogram, a vector of the same length as x.
#
smoother<-function(sm,x)
{
  lsm<-length(sm)
  lx<-length(x)
  
  #some error checking
  if (lsm>=lx)
  {
    stop("Error in smoother: length of smoother must be less than length of vector to be smoothed")
  }
  if (lsm%%2==0 || lsm==1)
  {
    stop("Error in smoother: sm must have odd length, at least 3")
  }
  
  #do the smoothing
  xaug<-c(NA,x,NA,x,NA,x)
  res<-numeric(lx)
  smbl<-(lsm-1)/2
  for (counter in 1:lx)
  {
    thisloc<-lx+2+counter
    h1<-xaug[(thisloc-smbl):(thisloc+smbl)]
    res[counter]<-sum(h1*sm,na.rm=TRUE)
  }
  
  return(res)
}

#The power spectrum, Brillinger's consistent estimator (5.6 of Brillinger's 2001 
#book).
#
#Args
#x                A time series as a vector
#detrend          TRUE or FALSE according to whether you want to linearly detrend 
#                   the time series first
#cutsym           TRUE or FALSE depending on whether you want to cut out the 
#                   symmetric part of the estimated spectrum in the output
#forvar           TRUE or FALSE depending on whether you want the spectrum rescaled
#                   so the mean (across all frequencies, including the symmetric part
#                   but not including the zero frequency) exactly equals the variance 
#                   of x (computed after detrending if that was done).
#BiasVariance     For adjusting the bias-variance tradeoff which comes from the 
#                   degree of smoothing selected
#
#Output - a list with these elements
#freq             The frequencies
#spec             The spectrum
#
#Details
#Frequency is here in units of cycles per sampling interval, whereas it was radians per 
#sampling interval in Brillinger.
#
myspecbrill<-function(x,detrend=TRUE,cutsym=TRUE,forvar=FALSE,BiasVariance=0.5)
{
  Tx<-length(x)
  
  #detrend, if desired
  if (detrend==TRUE)
  {
    tforx<-1:Tx
    x<-stats::residuals(stats::lm(x~tforx))
  } 
  
  #get the raw periodogram, with the zero frequency removed
  rawper<-myspecraw(x,detrend=FALSE,removezero=TRUE,cutsym=FALSE)
  I<-rawper$spec
  freq<-rawper$freq

  #now do the smoothing that makes the estimator
  BTx<-BiasVariance/sqrt(Tx) #adjust BT to adjust the bias-variance tradeoff
  TxBTx<-Tx*BTx
  fTxBTx<-floor(TxBTx)
  xforW<-(-fTxBTx:fTxBTx)/TxBTx
  WT<-(15/(16*2*pi))*((xforW-1)^2)*((xforW+1)^2)
  spec<-smoother(WT,I)
  spec<-2*pi*spec/(sqrt(Tx)*BiasVariance) #normalization

  #now rescale to make the mean exactly equal the variance of the original time series, if desired
  if (forvar==TRUE)
  {
    spec<-var(x)*spec/mean(spec)
  }
  
  #cut the symmetric part of the spectrum, if desired
  if (cutsym==TRUE)
  {
    spec<-spec[freq<=0.5]
    freq<-freq[freq<=0.5]
  }
  
  return(list(freq=freq,spec=spec))
}

#The above function is my attempt, begun 2020 11 19, to refactor the below function, largely inspired by 
#the need to make sure spectra appropriately integrate to equal variances even when smoothing is 
#used. But also I just saw some opportunities to clean things up. I kept the old function for testing.

#The power spectrum, Brillinger's consistent estimator (5.6 of Brillinger's 2001 
#book). The only difference from what is described there is frequency is here in 
#units of cycles per sampling interval in the output here, and was in radians per 
#sampling interval in Brillinger. Detrending is optionally done by this function. 
#The function actually returns the log scale power spectrum.
#
#Args
#x      A time series as a vector
#BiasVariance     For adjusting the bias-variance tradeoff which comes from the 
#                   degree of smoothing selected
#
myspecbrill_old<-function(x,detrend=TRUE,BiasVariance=0.5)
{
  Tx<-length(x)
  
  #detrend, if desired
  if (detrend==TRUE)
  {
    tforx<-1:Tx
    x<-residuals(lm(x~tforx))
  }
  
  #get the raw periodogram
  fftx<-fft(x)
  I<-(Mod(fftx))^2/(2*pi*Tx)
  I[1]<-0 #Set zero frequency to 0. Should be zero anyway, to within rounding error, if there was detrending
  freq<-(0:(Tx-1))/Tx
  freq<-2*pi*freq #to make frequencies be in units of radians per sampling interval
  
  #now do the smoothing that makes the estimator
  BTx<-BiasVariance/sqrt(Tx) #adjust BT to adjust the bias-variance tradeoff
  TxBTx<-Tx*BTx
  xforW<-(0:floor(TxBTx))/TxBTx
  WT<-(15/(16*2*pi))*((xforW-1)^2)*((xforW+1)^2)
  intWsquare<-(15/(32*pi))^2*4*pi*(1/9-4/7+6/5-4/3+1)
  spec<-WT[1]*I
  lenI<-length(I)
  for (AbsShift in 1:(length(WT)-1))
  {
    temp<-WT[AbsShift+1]*I
    #spec<-spec + temp([(AbsShift+1):end 1:AbsShift],:) + temp([(end-AbsShift+1):end 1:(end-AbsShift)],:);
    spec<-spec+temp[c((AbsShift+1):lenI,1:AbsShift)]+temp[c((lenI-AbsShift+1):lenI,1:(lenI-AbsShift))]
  }
  
  #remove the 0 frequency, and change the units back to cycles per sampling interval,
  #and cut the redundant part of the spectrum
  freq<-freq[-1]
  spec<-spec[-1]
  freq<-freq/(2*pi)
  spec<-spec[freq<=0.5]
  freq<-freq[freq<=0.5]
  
  #put in some normalization factors
  spec<-spec*2*pi/TxBTx
  
  #now get confidence intervals
  p<-0.95 #for 95% confidence intervals
  conf<-qnorm((1+p)/2,mean=0,sd=1)*0.4343*sqrt(2*pi*intWsquare/(TxBTx));
  
  return(list(freq=freq,log10spec=log10(spec),conf=conf))
}

#The spectral matrix, Brillinger's consistent estimator (7.4 of Brillinger's 2001 
#book).
#
#Args
#x                A vector-valued time series as a matrix, components by time (so time runs 
#                   along the rows). 
#detrend          TRUE or FALSE according to whether you want to linearly detrend 
#                   the time series first, each component time series separately.
#cutsym           TRUE or FALSE depending on whether you want to cut out the 
#                   frequencies (and the value that go with them) for the symmetric 
#                   part of the spectrum.
#forvar           TRUE or FALSE depending on whether you want things rescaled so the 
#                   means (across all frequencies, including the symmetric part but 
#                   not including the zero frequency) of the estimated spectra equal 
#                   the variance of the component time series (computed after detrending 
#                   if that was done); and the means of the estimated cospectra equal 
#                   the covariance of the component time series.
#BiasVariance     For adjusting the bias-variance tradeoff which comes from the 
#                   degree of smoothing selected
#
myspecmatbrill<-function(x,detrend=TRUE,cutsym=TRUE,forvar=FALSE,BiasVariance=0.5)
{
  Tx<-dim(x)[2] #length of time series
  N<-dim(x)[1] #number of components
  
  #detrend
  if (detrend==TRUE)
  {
    tforx<-1:Tx
    for (counter in 1:N)
    {
      y<-x[counter,]
      x[counter,]<-stats::residuals(stats::lm(y~tforx))
    }
  }
  
  rawspecmat<-myspecmatraw(x,detrend=FALSE,removezero=TRUE,cutsym=FALSE)
  freq<-rawspecmat$freq
  I<-rawspecmat$spec
  
  #now do the smoothing that makes the estimator
  BTx<-BiasVariance/sqrt(Tx) #adjust BT to adjust the bias-variance tradeoff
  TxBTx<-Tx*BTx
  fTxBTx<-floor(TxBTx)
  xforW<-(-fTxBTx:fTxBTx)/TxBTx
  WT<-(15/(16*2*pi))*((xforW-1)^2)*((xforW+1)^2)
  spec<-array(complex(real=0,imaginary=0),dim(I))
  for (a in 1:N)
  {
    for (b in 1:N)
    {
      spec[a,b,]<-smoother(WT,I[a,b,])
    }
  }
  spec<-2*pi*spec/(sqrt(Tx)*BiasVariance) #normalization
  
  #now rescale to make the means equal the covariances of the original time series, if desired
  if (forvar==TRUE)
  {
    for (a in 1:N)
    {
      for (b in 1:N)
      {
        spec[a,b,]<-cov(x[a,],x[b,])*spec[a,b,]/mean(spec[a,b,])
      }
    }
  }
  
  #cut the symmetric part of the spectrum, if desired
  if (cutsym==TRUE)
  {
    spec<-spec[,,freq<=0.5]
    freq<-freq[freq<=0.5]
  }
  
  return(list(freq=freq,spec=spec))
}

#The above function is my attempt, begun 2020 11 19, to refactor the below function, largely inspired by 
#the need to make sure spectra appropriately integrate to equal variances even when smoothing is 
#used. But also I just saw some opportunities to clean things up. I kept the old function for testing.

#The spectral matrix, Brillinger's consistent estimator (7.4 of Brillinger's 2001 
#book). The only difference from what is described there is frequency is here in 
#units of cycles per sampling interval in the output, and was in radians per sampling 
#interval in Brillinger. Linear detrending (and de-meaning) of the individual 
#component time series is optionally done.
#
#Args
#x      A vector-valued time series as a matrix, components by time (so time runs 
#         along the rows). 
#BiasVariance     For adjusting the bias-variance tradeoff which comes from the 
#                   degree of smoothing selected
#
myspecmatbrill_old<-function(x,detrend=TRUE,BiasVariance=0.5)
{
  Tx<-dim(x)[2] #length of time series
  N<-dim(x)[1] #number of components
  
  #detrend
  if (detrend)
  {
    tforx<-1:Tx
    for (counter in 1:N)
    {
      y<-x[counter,]
      x[counter,]<-residuals(lm(y~tforx))  
    }
  }
  
  #get the raw periodogram
  fftx<-mymvfft(x)
  I<-array(complex(real=0,imaginary=0),dim=c(N,N,Tx))
  for (a in 1:N)
  {
    for (b in 1:N)
    {
      I[a,b,]<-fftx[a,]*Conj(fftx[b,])
    }
  }
  I<-I/(2*pi*Tx)
  freq<-(0:(Tx-1))/Tx
  freq<-2*pi*freq #to make frequencies be in units of radians per sampling interval, for now
  
  #now do the smoothing that makes the estimator
  BTx<-BiasVariance/sqrt(Tx) #adjust BT to adjust the bias-variance tradeoff
  TxBTx<-Tx*BTx
  xforW<-(0:floor(TxBTx))/TxBTx
  WT<-(15/(16*2*pi))*((xforW-1)^2)*((xforW+1)^2)
  intWsquare<-(15/(32*pi))^2*4*pi*(1/9-4/7+6/5-4/3+1)
  specmat<-WT[1]*I
  lenI<-dim(I)[3]
  for (AbsShift in 1:(length(WT)-1))
  {
    temp<-WT[AbsShift+1]*I
    specmat<-specmat+temp[,,c((AbsShift+1):lenI,1:AbsShift)]+temp[,,c((lenI-AbsShift+1):lenI,1:(lenI-AbsShift))]
  }
  
  #remove the 0 frequency, and change the units back to radians per sampling interval,
  #and cut the redundant part 
  freq<-freq[-1]
  specmat<-specmat[,,-1,drop=FALSE]
  freq<-freq/(2*pi)
  specmat<-specmat[,,freq<=0.5]
  freq<-freq[freq<=0.5]
  
  #put in some normalization factors
  specmat<-specmat*2*pi/TxBTx
  
  return(list(freq=freq,specmat=specmat))
}
