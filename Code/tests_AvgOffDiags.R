context("AvgOffDiags.R")

source("AvgOffDiags.R")

test_that("test it",{
  #***for the default case
  
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
  
  #***for the case where resc is provided
  
  #test error catching
  resc<-c(1,2,3)
  testthat::expect_error(avgoffdiags(A,resc),"Error in avgoffdiags: bad value for resc")
  resc<-c(1,2,3,4,NA)
  testthat::expect_error(avgoffdiags(A,resc),"Error in avgoffdiags: bad value for resc")
  
  #test a case that should equal the default case
  resc<-rep(1,5)
  testthat::expect_equal(avgoffdiags(A),avgoffdiags(A,resc))
  
  #now test a case which actually uses a non-trivial resc
  resc<-1:5
  res1<-avgoffdiags(A,resc)
  rescmat<-outer(sqrt(resc),sqrt(resc),FUN="*")
  rescmat<-array(rescmat,c(5,5,3))
  A<-A/rescmat
  for (counter in 1:3)
  {
    h<-A[,,counter]
    diag(h)<-NA
    A[,,counter]<-h
  }
  res2<-apply(FUN=mean,X=A,MARGIN=3,na.rm=TRUE)
  testthat::expect_equal(res1,res2)
})