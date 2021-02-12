context("AvgOffDiags.R")

source("AvgOffDiags.R")

test_that("test it",{
  #test error catching
  A<-array(1,c(2,3,4))
  testthat::expect_error(avgoffdiags(A),"Error in avgoffdiags: A must have dimensions n by n by m")
  A<-array(1,c(2,2,2))
  A[1,1,1]<-NA
  testthat::expect_error(avgoffdiags(A),"Error in avgoffdiags: A must have all finite entries")

  #test some cases
  set.seed(101)
  A<-array(rnorm(75),c(5,5,3))
  res<-avgoffdiags(A) 
  
  h<-A[,,1]
  diag(h)<-NA
  testthat::expect_equal(res[1],mean(h,na.rm=TRUE))

  h<-A[,,2]
  diag(h)<-NA
  testthat::expect_equal(res[2],mean(h,na.rm=TRUE))

  h<-A[,,3]
  diag(h)<-NA
  testthat::expect_equal(res[3],mean(h,na.rm=TRUE))
})