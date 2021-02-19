context("MnCoefDet_Tiny2.R")

source("MnCoefDet_Tiny2.R")

test_that("test it",{
  #test error catching
  testthat::expect_error(mncoefdet_tiny2(kelp=array(1,c(2,2)),NO3=array(2,c(2,3)),waves=array(3,c(2,2))),
                         "Error in mncoefdet_tiny2: kelp, NO3 and waves input matrices must be the same dimensions")
  testthat::expect_error(mncoefdet_tiny2(kelp=array(1,c(2,2)),NO3=array(2,c(2,2)),waves=array(3,c(2,3))),
                         "Error in mncoefdet_tiny2: kelp, NO3 and waves input matrices must be the same dimensions")

  #now test some cases with legit inputs
  set.seed(101)
  kelp<-matrix(rnorm(1000),10,100)
  NO3<-matrix(rnorm(1000),10,100)
  waves<-matrix(rnorm(1000),10,100)
  
  res<-mncoefdet_tiny2(kelp,NO3,waves)
  res2<-numeric(5)
  for (counter in 1:10)
  {
    m<-lm(kelp[counter,2:100]~NO3[counter,2:100]+NO3[counter,1:99]+waves[counter,2:100]+waves[counter,1:99])
    res2[1]<-res2[1]+summary(m)$r.squared
    res2[2:3]<-res2[2:3]+unname(coef(m)[2:3])
    res2[4:5]<-res2[4:5]+unname(coef(m)[4:5])
  }
  testthat::expect_equal(res,res2/10)
})