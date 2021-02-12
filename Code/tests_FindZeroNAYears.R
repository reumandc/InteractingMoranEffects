context("FillZeroNAYears.R")

source("FindZeroNAYears.R")

test_that("test it",{
  #error catching
  x<-c(1,2,3,4,5,6,7)
  testthat::expect_error(find_zero_NA_years(x),"Error in find_zero_NA_years: x needs to have length divisible by four")  
  
  #test some cases
  x<-c(1,2,3,4,1,2,3,4)
  res<-find_zero_NA_years(x)
  testthat::expect_equal(res,0)
  
  x<-c(1,2,3,4,1,2,3,4,0,0,0,0)
  res<-find_zero_NA_years(x)
  testthat::expect_equal(res,1)

  x<-c(1,2,3,4,NA,0,1,2,1,2,3,4)
  res<-find_zero_NA_years(x)
  testthat::expect_equal(res,0)

  x<-c(NA,1,2,3,1,2,3,4,1,2,3,4)
  res<-find_zero_NA_years(x)
  testthat::expect_equal(res,0)

  x<-c(1,2,3,4,1,2,3,4,0,0,0,1)
  res<-find_zero_NA_years(x)
  testthat::expect_equal(res,0)

  x<-c(1,2,3,4,1,2,3,4,0,0,0,0)
  res<-find_zero_NA_years(x)
  testthat::expect_equal(res,1)

  x<-c(1,2,3,4,1,2,0,NA,0,NA,1,2,1,2,3,4)
  res<-find_zero_NA_years(x)
  testthat::expect_equal(res,0)
})