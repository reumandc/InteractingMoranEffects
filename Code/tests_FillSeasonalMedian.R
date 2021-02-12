context("FillSeasonalMedian.R")

source("../Code/FillSeasonalMedian.R")

test_that("test it",{
  set.seed(101)
  x<-rnorm(100)
  x[1]<-NA
  x[6]<-NA
  x[11]<-NA
  x[16]<-NA
  res<-fill_with_seasonal_median(x)
  expect_equal(x[-c(1,6,11,16)],res[-c(1,6,11,16)])
  expect_equal(res[1],median(x[1:length(x) %% 4==1],na.rm=TRUE))
  expect_equal(res[6],median(x[1:length(x) %% 4==2],na.rm=TRUE))
  expect_equal(res[11],median(x[1:length(x) %% 4==3],na.rm=TRUE))
  expect_equal(res[16],median(x[1:length(x) %% 4==0],na.rm=TRUE))
})


