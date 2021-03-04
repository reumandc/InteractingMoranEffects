context("PlotPhaseFunc.R")

source("PlotPhaseFunc.R")

test_that("test it",{
  #this just tests the error handling, tests of the plotting capabilities are in plottests_PlotPhaseFunc.R
  testthat::expect_error(plotphasefunc(1:4,1:4),"Error in plotphasefunc: z must be a complex vector")
  testthat::expect_error(plotphasefunc("test",complex(real=1:4,imaginary=1:4)),"Error in plotphasefunc: t must be a numeric vector")  
  testthat::expect_error(plotphasefunc(1:3,complex(real=1:4,imaginary=1:4)),"Error in plotphasefunc: t and z arguments must be vectors of the same length")  
  testthat::expect_error(plotphasefunc(3:1,complex(real=1:3,imaginary=1:3)),"Error in plotphasefunc: values of t must be strictly increasing")
})