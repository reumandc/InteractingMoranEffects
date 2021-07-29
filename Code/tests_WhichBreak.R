context("WhichBreak.R")

source("WhichBreak.R")

test_that("test it",{
  #test the error catching
  testthat::expect_error(WhichBreak(c(1,2,4),rep(.5,5)),"Error in WhichBreak: breaks must have evenly spaced values")
  testthat::expect_error(WhichBreak(c(3,2,1),rep(.5,5)),"Error in WhichBreak: breaks must have increasing values")
  testthat::expect_error(WhichBreak(c(1,2,3),c(-1,2,2.5)),"Error in WhichBreak: some entries of vals fall outside the range of breaks")
  
  #test a simple case
  testthat::expect_equal(WhichBreak(c(0,1,2,3,4),c(.1,1.1,2,3,4)),c(1,2,2,3,4))
})