context("MnCor_Tiny.R")

source("MnCor_Tiny.R")

test_that("test it",{
  #test error catching
  m1<-matrix(rnorm(10),2,5)
  m2<-matrix(rnorm(10),5,2)
  testthat::expect_error(mncor_tiny(m1,m2,0),"Error in mncor_tiny: m1 and m2 must have the same dimension")
  m2<-t(m2)
  testthat::expect_error(mncor_tiny(m1,m2,-1),"Error in mncor_tiny: lag must be >=0")
  
  #test some cases
  set.seed(101)
  m1<-matrix(rnorm(200),c(2,100))
  m2<-array(rnorm(200),c(2,100))
  
  res<-mncor_tiny(m1,m2,0)
  res2<-(cor(m1[1,],m2[1,])+cor(m1[2,],m2[2,]))/2
  testthat::expect_equal(res,res2)
  
  res<-mncor_tiny(m1,m2,1)
  res2<-(cor(m1[1,2:100],m2[1,1:99])+cor(m1[2,2:100],m2[2,1:99]))/2
  testthat::expect_equal(res,res2)
})