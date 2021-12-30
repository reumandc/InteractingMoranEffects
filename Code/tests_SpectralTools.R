context("SpectralTools.R")

source("../Code/SpectralTools.R")

test_that("test the raw periodogram function, myspecraw",{
  set.seed(101)
  
  x<-rnorm(100)
  resTTT<-myspecraw(x,detrend=TRUE,removezero=TRUE,cutsym=TRUE)
  resTTF<-myspecraw(x,detrend=TRUE,removezero=TRUE,cutsym=FALSE)
  resTFT<-myspecraw(x,detrend=TRUE,removezero=FALSE,cutsym=TRUE)
  resTFF<-myspecraw(x,detrend=TRUE,removezero=FALSE,cutsym=FALSE)
  resFTT<-myspecraw(x,detrend=FALSE,removezero=TRUE,cutsym=TRUE)
  resFTF<-myspecraw(x,detrend=FALSE,removezero=TRUE,cutsym=FALSE)
  resFFT<-myspecraw(x,detrend=FALSE,removezero=FALSE,cutsym=TRUE)
  resFFF<-myspecraw(x,detrend=FALSE,removezero=FALSE,cutsym=FALSE)
  
  #check output is the right format
  expect_equal(class(resTTT),"list")
  expect_equal(class(resTTF),"list")
  expect_equal(class(resTFT),"list")
  expect_equal(class(resTFF),"list")
  expect_equal(class(resFTT),"list")
  expect_equal(class(resFTF),"list")
  expect_equal(class(resFFT),"list")
  expect_equal(class(resFFF),"list")
  
  expect_equal(names(resTTT),c("freq","spec"))
  expect_equal(names(resTTF),c("freq","spec"))
  expect_equal(names(resTFT),c("freq","spec"))
  expect_equal(names(resTFF),c("freq","spec"))
  expect_equal(names(resFTT),c("freq","spec"))
  expect_equal(names(resFTF),c("freq","spec"))
  expect_equal(names(resFFT),c("freq","spec"))
  expect_equal(names(resFFF),c("freq","spec"))
  
  frFF<-(0:(length(x)-1))/length(x)
  frFT<-frFF[frFF<=0.5]
  frTF<-(1:(length(x)-1))/length(x)
  frTT<-frTF[frTF<=0.5]
  expect_equal(resTTT$freq,frTT)
  expect_equal(resTTF$freq,frTF)
  expect_equal(resTFT$freq,frFT)
  expect_equal(resTFF$freq,frFF)
  expect_equal(resFTT$freq,frTT)
  expect_equal(resFTF$freq,frTF)
  expect_equal(resFFT$freq,frFT)
  expect_equal(resFFF$freq,frFF)
  
  expect_equal(length(resTTT$freq),length(resTTT$spec))
  expect_equal(length(resTTF$freq),length(resTTF$spec))
  expect_equal(length(resTFT$freq),length(resTFT$spec))
  expect_equal(length(resTFF$freq),length(resTFF$spec))
  expect_equal(length(resFTT$freq),length(resFTT$spec))
  expect_equal(length(resFTF$freq),length(resFTF$spec))
  expect_equal(length(resFFT$freq),length(resFFT$spec))
  expect_equal(length(resFFF$freq),length(resFFF$spec))
  
  #now check the average across the spectrum is the variance
  expect_true(abs(var(x)-mean(resFTF$spec))<1e-15)

  #same thing again for a different length time series
  x<-rnorm(101)
  resTTT<-myspecraw(x,detrend=TRUE,removezero=TRUE,cutsym=TRUE)
  resTTF<-myspecraw(x,detrend=TRUE,removezero=TRUE,cutsym=FALSE)
  resTFT<-myspecraw(x,detrend=TRUE,removezero=FALSE,cutsym=TRUE)
  resTFF<-myspecraw(x,detrend=TRUE,removezero=FALSE,cutsym=FALSE)
  resFTT<-myspecraw(x,detrend=FALSE,removezero=TRUE,cutsym=TRUE)
  resFTF<-myspecraw(x,detrend=FALSE,removezero=TRUE,cutsym=FALSE)
  resFFT<-myspecraw(x,detrend=FALSE,removezero=FALSE,cutsym=TRUE)
  resFFF<-myspecraw(x,detrend=FALSE,removezero=FALSE,cutsym=FALSE)
  
  #check output is the right format
  expect_equal(class(resTTT),"list")
  expect_equal(class(resTTF),"list")
  expect_equal(class(resTFT),"list")
  expect_equal(class(resTFF),"list")
  expect_equal(class(resFTT),"list")
  expect_equal(class(resFTF),"list")
  expect_equal(class(resFFT),"list")
  expect_equal(class(resFFF),"list")
  
  expect_equal(names(resTTT),c("freq","spec"))
  expect_equal(names(resTTF),c("freq","spec"))
  expect_equal(names(resTFT),c("freq","spec"))
  expect_equal(names(resTFF),c("freq","spec"))
  expect_equal(names(resFTT),c("freq","spec"))
  expect_equal(names(resFTF),c("freq","spec"))
  expect_equal(names(resFFT),c("freq","spec"))
  expect_equal(names(resFFF),c("freq","spec"))
  
  frFF<-(0:(length(x)-1))/length(x)
  frFT<-frFF[frFF<=0.5]
  frTF<-(1:(length(x)-1))/length(x)
  frTT<-frTF[frTF<=0.5]
  expect_equal(resTTT$freq,frTT)
  expect_equal(resTTF$freq,frTF)
  expect_equal(resTFT$freq,frFT)
  expect_equal(resTFF$freq,frFF)
  expect_equal(resFTT$freq,frTT)
  expect_equal(resFTF$freq,frTF)
  expect_equal(resFFT$freq,frFT)
  expect_equal(resFFF$freq,frFF)
  
  expect_equal(length(resTTT$freq),length(resTTT$spec))
  expect_equal(length(resTTF$freq),length(resTTF$spec))
  expect_equal(length(resTFT$freq),length(resTFT$spec))
  expect_equal(length(resTFF$freq),length(resTFF$spec))
  expect_equal(length(resFTT$freq),length(resFTT$spec))
  expect_equal(length(resFTF$freq),length(resFTF$spec))
  expect_equal(length(resFFT$freq),length(resFFT$spec))
  expect_equal(length(resFFF$freq),length(resFFF$spec))
  
  #now check the average across the spectrum is the variance
  expect_true(abs(var(x)-mean(resFTF$spec))<1e-15)
})

test_that("test myspecmatraw",{
  #test the format of the output
  set.seed(101)
  x<-matrix(rnorm(1000),10,100)
  resTTT<-myspecmatraw(x,detrend=TRUE,removezero=TRUE,cutsym=TRUE)
  resTTF<-myspecmatraw(x,detrend=TRUE,removezero=TRUE,cutsym=FALSE)
  resTFT<-myspecmatraw(x,detrend=TRUE,removezero=FALSE,cutsym=TRUE)
  resTFF<-myspecmatraw(x,detrend=TRUE,removezero=FALSE,cutsym=FALSE)
  resFTT<-myspecmatraw(x,detrend=FALSE,removezero=TRUE,cutsym=TRUE)
  resFTF<-myspecmatraw(x,detrend=FALSE,removezero=TRUE,cutsym=FALSE)
  resFFT<-myspecmatraw(x,detrend=FALSE,removezero=FALSE,cutsym=TRUE)
  resFFF<-myspecmatraw(x,detrend=FALSE,removezero=FALSE,cutsym=FALSE)
  
  expect_equal(class(resTTT),"list")
  expect_equal(class(resTFT),"list")
  expect_equal(class(resFTT),"list")
  expect_equal(class(resFFT),"list")
  expect_equal(class(resTTF),"list")
  expect_equal(class(resTFF),"list")
  expect_equal(class(resFTF),"list")
  expect_equal(class(resFFF),"list")
  
  expect_equal(names(resTTT),c("freq","spec"))
  expect_equal(names(resTFT),c("freq","spec"))
  expect_equal(names(resFTT),c("freq","spec"))
  expect_equal(names(resFFT),c("freq","spec"))
  expect_equal(names(resTTF),c("freq","spec"))
  expect_equal(names(resTFF),c("freq","spec"))
  expect_equal(names(resFTF),c("freq","spec"))
  expect_equal(names(resFFF),c("freq","spec"))
  
  lx<-dim(x)[2]
  frFF<-(0:(lx-1))/lx
  frTF<-(1:(lx-1))/lx
  frFT<-frFF[frFF<=0.5]
  frTT<-frTF[frTF<=0.5]
  expect_equal(resTTT$freq,frTT)
  expect_equal(resTFT$freq,frFT)
  expect_equal(resFTT$freq,frTT)
  expect_equal(resFFT$freq,frFT)
  expect_equal(resTTF$freq,frTF)
  expect_equal(resTFF$freq,frFF)
  expect_equal(resFTF$freq,frTF)
  expect_equal(resFFF$freq,frFF)
  
  expect_equal(length(resTTT$freq),dim(resTTT$spec)[3])  
  expect_equal(length(resTFT$freq),dim(resTFT$spec)[3])
  expect_equal(length(resFTT$freq),dim(resFTT$spec)[3])
  expect_equal(length(resFFT$freq),dim(resFFT$spec)[3])
  expect_equal(length(resTTF$freq),dim(resTTF$spec)[3])  
  expect_equal(length(resTFF$freq),dim(resTFF$spec)[3])
  expect_equal(length(resFTF$freq),dim(resFTF$spec)[3])
  expect_equal(length(resFFF$freq),dim(resFFF$spec)[3])
  
  #now check the mean across the cospectrum is the covariance
  N<-dim(x)[1]
  for (a in 1:N)
  {
    for (b in 1:N)
    {
      expect_true(abs(cov(x[a,],x[b,])-mean(Re(resFTF$spec[a,b,])))<1e-15)
    }
  }
  xdetr<-x
  tforx<-1:lx
  for (counter in 1:N)
  {
    y<-xdetr[counter,]
    xdetr[counter,]<-residuals(lm(y~tforx))  
  }
  for (a in 1:N)
  {
    for (b in 1:N)
    {
      expect_true(abs(cov(xdetr[a,],xdetr[b,])-mean(Re(resTTF$spec[a,b,])))<1e-15)
    }
  }
})

test_that("test the mymvfft and imymvfft convenience functions",{
  set.seed(101)
  
  x<-matrix(rnorm(1000),10,100)
  res<-mymvfft(x)
  ix<-imymvfft(res)
  
  expect_equal(res[1,],fft(x[1,]))
  expect_equal(res[2,],fft(x[2,]))
  expect_equal(res[3,],fft(x[3,]))

  expect_equal(Im(ix),matrix(0,nrow(x),ncol(x)))  
  expect_equal(x,Re(ix))
})

test_that("test the smoother function",{
  #test error catching
  expect_error(smoother(1:10,1:10),"Error in smoother: length of smoother must be less than length of vector to be smoothed")
  expect_error(smoother(1,1:10),"Error in smoother: sm must have odd length, at least 3")
  expect_error(smoother(1:4,1:10),"Error in smoother: sm must have odd length, at least 3")

  #test on an example  
  x<-1:10
  sm<-c(.4,1,.5)
  res<-smoother(sm,x)
  xaug<-c(0,x)
  manres<-xaug*sm[2]+xaug[c(11,1:10)]*sm[1]+xaug[c(2:11,1)]*sm[3]
  expect_equal(res,manres[2:11])
  
  #test on another example with a longer smoother and a longer vector
  set.seed(101)
  x<-rnorm(100)
  sm<-c(.1,.2,.3,.2,.1)
  res<-smoother(sm,x)
  xaug<-c(0,x)
  manres<-xaug*sm[3]+
    xaug[c(100:101,1:99)]*sm[1]+
    xaug[c(101,1:100)]*sm[2]+
    xaug[c(2:101,1)]*sm[4]+
    xaug[c(3:101,1:2)]*sm[5]
  expect_equal(res,manres[2:101])
})

test_that("test myspecbrill",{
  #test the format of the output
  set.seed(101)
  x<-rnorm(100)
  resTTT<-myspecbrill(x,detrend=TRUE,cutsym=TRUE,forvar=TRUE,BiasVariance=0.5)
  resTFT<-myspecbrill(x,detrend=TRUE,cutsym=FALSE,forvar=TRUE,BiasVariance=0.5)
  resFTT<-myspecbrill(x,detrend=FALSE,cutsym=TRUE,forvar=TRUE,BiasVariance=0.5)
  resFFT<-myspecbrill(x,detrend=FALSE,cutsym=FALSE,forvar=TRUE,BiasVariance=0.5)
  resTTF<-myspecbrill(x,detrend=TRUE,cutsym=TRUE,forvar=FALSE,BiasVariance=0.5)
  resTFF<-myspecbrill(x,detrend=TRUE,cutsym=FALSE,forvar=FALSE,BiasVariance=0.5)
  resFTF<-myspecbrill(x,detrend=FALSE,cutsym=TRUE,forvar=FALSE,BiasVariance=0.5)
  resFFF<-myspecbrill(x,detrend=FALSE,cutsym=FALSE,forvar=FALSE,BiasVariance=0.5)

  expect_equal(class(resTTT),"list")
  expect_equal(class(resTFT),"list")
  expect_equal(class(resFTT),"list")
  expect_equal(class(resFFT),"list")
  expect_equal(class(resTTF),"list")
  expect_equal(class(resTFF),"list")
  expect_equal(class(resFTF),"list")
  expect_equal(class(resFFF),"list")
  
  expect_equal(names(resTTT),c("freq","spec"))
  expect_equal(names(resTFT),c("freq","spec"))
  expect_equal(names(resFTT),c("freq","spec"))
  expect_equal(names(resFFT),c("freq","spec"))
  expect_equal(names(resTTF),c("freq","spec"))
  expect_equal(names(resTFF),c("freq","spec"))
  expect_equal(names(resFTF),c("freq","spec"))
  expect_equal(names(resFFF),c("freq","spec"))
  
  frF<-(1:(length(x)-1))/length(x)
  frT<-frF[frF<=0.5]
  expect_equal(resTTT$freq,frT)
  expect_equal(resTFT$freq,frF)
  expect_equal(resFTT$freq,frT)
  expect_equal(resFFT$freq,frF)
  expect_equal(resTTF$freq,frT)
  expect_equal(resTFF$freq,frF)
  expect_equal(resFTF$freq,frT)
  expect_equal(resFFF$freq,frF)
  
  expect_equal(length(resTTT$freq),length(resTTT$spec))  
  expect_equal(length(resTFT$freq),length(resTFT$spec))
  expect_equal(length(resFTT$freq),length(resFTT$spec))
  expect_equal(length(resFFT$freq),length(resFFT$spec))
  expect_equal(length(resTTF$freq),length(resTTF$spec))  
  expect_equal(length(resTFF$freq),length(resTFF$spec))
  expect_equal(length(resFTF$freq),length(resFTF$spec))
  expect_equal(length(resFFF$freq),length(resFFF$spec))
  
  #comparison with the old version of the code
  res_old<-myspecbrill_old(x,detrend=TRUE,BiasVariance=0.5)
  res_new<-resTTF
  expect_equal(res_old$freq,res_new$freq)
  expect_true(max(abs((10^res_old$log10spec)*2*pi-res_new$spec))<1e-15)
  
  #now check the average across the spectrum is the variance, exactly, when forvar==TRUE
  expect_true(abs(var(x)-mean(resFFT$spec))<1e-15)
  tforx<-1:length(x)
  xdetr<-stats::residuals(stats::lm(x~tforx))
  expect_true(abs(var(xdetr)-mean(resTFT$spec))<1e-15)
  
  #check the average across the spectrum is the variance, to within, say, 5%, even 
  #when forvar==FALSE
  expect_true(abs((var(x)-mean(resFFF$spec))/var(x))<.05)
})

test_that("test myspecmatbrill",{
  #test the format of the output
  set.seed(101)
  sig<-matrix(.3,10,10)
  diag(sig)<-1
  x<-t(mvtnorm::rmvnorm(100,mean=rep(0,10),sigma=sig))

  resTTT<-myspecmatbrill(x,detrend=TRUE,cutsym=TRUE,forvar=TRUE,BiasVariance=0.5)
  resTTF<-myspecmatbrill(x,detrend=TRUE,cutsym=TRUE,forvar=FALSE,BiasVariance=0.5)
  resTFT<-myspecmatbrill(x,detrend=TRUE,cutsym=FALSE,forvar=TRUE,BiasVariance=0.5)
  resTFF<-myspecmatbrill(x,detrend=TRUE,cutsym=FALSE,forvar=FALSE,BiasVariance=0.5)
  resFTT<-myspecmatbrill(x,detrend=FALSE,cutsym=TRUE,forvar=TRUE,BiasVariance=0.5)
  resFTF<-myspecmatbrill(x,detrend=FALSE,cutsym=TRUE,forvar=FALSE,BiasVariance=0.5)
  resFFT<-myspecmatbrill(x,detrend=FALSE,cutsym=FALSE,forvar=TRUE,BiasVariance=0.5)
  resFFF<-myspecmatbrill(x,detrend=FALSE,cutsym=FALSE,forvar=FALSE,BiasVariance=0.5)

  expect_equal(class(resTTT),"list")
  expect_equal(class(resTFT),"list")
  expect_equal(class(resFTT),"list")
  expect_equal(class(resFFT),"list")
  expect_equal(class(resTTF),"list")
  expect_equal(class(resTFF),"list")
  expect_equal(class(resFTF),"list")
  expect_equal(class(resFFF),"list")
  
  expect_equal(names(resTTT),c("freq","spec"))
  expect_equal(names(resTFT),c("freq","spec"))
  expect_equal(names(resFTT),c("freq","spec"))
  expect_equal(names(resFFT),c("freq","spec"))
  expect_equal(names(resTTF),c("freq","spec"))
  expect_equal(names(resTFF),c("freq","spec"))
  expect_equal(names(resFTF),c("freq","spec"))
  expect_equal(names(resFFF),c("freq","spec"))
  
  lx<-dim(x)[2]
  frF<-(1:(lx-1))/lx
  frT<-frF[frF<=0.5]
  expect_equal(resTTT$freq,frT)
  expect_equal(resTFT$freq,frF)
  expect_equal(resFTT$freq,frT)
  expect_equal(resFFT$freq,frF)
  expect_equal(resTTF$freq,frT)
  expect_equal(resTFF$freq,frF)
  expect_equal(resFTF$freq,frT)
  expect_equal(resFFF$freq,frF)
  
  expect_equal(length(resTTT$freq),dim(resTTT$spec)[3])  
  expect_equal(length(resTFT$freq),dim(resTFT$spec)[3])
  expect_equal(length(resFTT$freq),dim(resFTT$spec)[3])
  expect_equal(length(resFFT$freq),dim(resFFT$spec)[3])
  expect_equal(length(resTTF$freq),dim(resTTF$spec)[3])  
  expect_equal(length(resTFF$freq),dim(resTFF$spec)[3])
  expect_equal(length(resFTF$freq),dim(resFTF$spec)[3])
  expect_equal(length(resFFF$freq),dim(resFFF$spec)[3])
  
  #comparison with the old version of the code
  res_old<-myspecmatbrill_old(x,detrend=TRUE,BiasVariance=0.5)
  res_new<-resTTF
  expect_equal(res_old$freq,res_new$freq)
  expect_true(max(abs(res_old$spec*2*pi-res_new$spec))<1e-15)
  
  #now check the average across the cospectrum is exactly the covariance, when forvar=TRUE
  N<-dim(x)[1]
  for (a in 1:N)
  {
    for (b in 1:N)
    {
      expect_true(abs(cov(x[a,],x[b,])-mean(resFFT$spec[a,b,]))<1e-15)
    }
  }
  xdetr<-x
  tforx<-1:lx
  for (counter in 1:N)
  {
    y<-xdetr[counter,]
    xdetr[counter,]<-residuals(lm(y~tforx))  
  }
  for (a in 1:N)
  {
    for (b in 1:N)
    {
      expect_true(abs(cov(xdetr[a,],xdetr[b,])-mean(resTFT$spec[a,b,]))<1e-15)
    }
  }
  
  #now check that averages across cospectrum are covariances to within 5% even when forvar=FALSE
  for (a in 1:N)
  {
    for (b in 1:N)
    {
      expect_true(abs((cov(x[a,],x[b,])-mean(resFFF$spec[a,b,]))/cov(x[a,],x[b,]))<.05)
    }
  }
  
  #check the same thing again with BiasVariance taking a different value
  resFFF<-myspecmatbrill(x,detrend=FALSE,cutsym=FALSE,forvar=FALSE,BiasVariance=1)
  for (a in 1:N)
  {
    for (b in 1:N)
    {
      expect_true(abs((cov(x[a,],x[b,])-mean(resFFF$spec[a,b,]))/cov(x[a,],x[b,]))<.05)
    }
  }
})
