context("IntroTheoryIA.R")

source("IntroTheoryIA.R")

test_that("test the error catching",{
  a<-1
  b<-1
  covmat<-diag(4)
  covmat[1,1]<-.9
  testthat::expect_error(IntroTheoryIA(a,b,covmat),"Error in IntroTheoryIA: covmat should have diagonal entries equal to 1")
})

test_that("test for the values given in Sheppard et al (2019), Plos Comp Biol, SI section S1",{
  a<-1
  b<-1

  covmat<-matrix(1,4,4)  
  covmat[1,4]<-1.1
  covmat[4,1]<-1.1
  h<-IntroTheoryIA(a,b,covmat)
  testthat::expect_true(is.na(h['sync']))
  
  covmat<-matrix(c(1,.9,.5,.5,
                   .9,1,.5,.5,
                   .5,.5,1,.5,
                   .5,.5,.5,1),4,4,byrow = TRUE)
  testthat::expect_equal(IntroTheoryIA(a,b,covmat),c(sync=0.8,sync.ni=0.7))
})
  