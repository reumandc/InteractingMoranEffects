context("PickLags.R")

source("PickLags.R")

test_that("test it",{
  
  #***test the function makedf
  
  #make up some data and call the function
  set.seed(101)
  T<-10
  y<-rnorm(T)
  x1<-rnorm(T)
  x2<-rnorm(T)
  maxlags<-c(3,2,1)
  dat<-makedf(y,x1,x2,maxlags)

  #test the column names and dimensions of the result
  testthat::expect_equal(names(dat),c("y","y_l1","y_l2","y_l3","x1_l0","x1_l1","x1_l2","x2_l0","x2_l1"))
  testthat::expect_equal(dim(dat)[1],7)

  #test the response variable column
  testthat::expect_equal(dat$y,y[4:10])

  #test the AR columns
  testthat::expect_equal(dat$y_l1,as.vector(y[3:9]))
  testthat::expect_equal(dat$y_l2,as.vector(y[2:8]))
  testthat::expect_equal(dat$y_l3,as.vector(y[1:7]))
  
  #test the x1 columns
  testthat::expect_equal(dat$x1_l0,x1[4:10])
  testthat::expect_equal(dat$x1_l1,x1[3:9])
  testthat::expect_equal(dat$x1_l2,x1[2:8])
  
  #test the x2 columns
  testthat::expect_equal(dat$x2_l0,x2[4:10])
  testthat::expect_equal(dat$x2_l1,x2[3:9])
  
  #***test the function makeformula
  
  #test a generic case
  res<-makeformula(3,2,1)
  testthat::expect_equal(class(res),"formula")
  testthat::expect_equal(res,as.formula("y~y_l1+y_l2+y_l3+x1_l0+x1_l1+x1_l2+x2_l0+x2_l1-1"))
  
  #test some boundary cases - one lag set at a null value
  res<-makeformula(0,2,1)
  testthat::expect_equal(class(res),"formula")
  testthat::expect_equal(res,as.formula("y~x1_l0+x1_l1+x1_l2+x2_l0+x2_l1-1"))

  res<-makeformula(2,-1,1)
  testthat::expect_equal(class(res),"formula")
  testthat::expect_equal(res,as.formula("y~y_l1+y_l2+x2_l0+x2_l1-1"))
  
  res<-makeformula(2,0,-1)
  testthat::expect_equal(class(res),"formula")
  testthat::expect_equal(res,as.formula("y~y_l1+y_l2+x1_l0-1"))
  
  #test some boundary cases - two lags set at null values
  res<-makeformula(0,-1,1)
  testthat::expect_equal(class(res),"formula")
  testthat::expect_equal(res,as.formula("y~x2_l0+x2_l1-1"))

  res<-makeformula(1,-1,-1)
  testthat::expect_equal(class(res),"formula")
  testthat::expect_equal(res,as.formula("y~y_l1-1"))
  
  res<-makeformula(0,1,-1)
  testthat::expect_equal(class(res),"formula")
  testthat::expect_equal(res,as.formula("y~x1_l0+x1_l1-1"))
  
  #***test the function picklags_vect
  
  #make up some data
  set.seed(101)
  T<-100
  y<-rnorm(T)
  x1<-rnorm(T)
  x2<-rnorm(T)
  maxlags<-c(3,2,1)

  #run the function
  res<-picklags_vect(y,x1,x2,maxlags)  

  #do some of the regressions "by hand" to check
  ind<-4
  dat<-makedf(y,x1,x2,maxlags)  
  formula<-makeformula(lagAR=res$lagAR[ind],lag1=res$lag1[ind],lag2=res$lag2[ind])
  mod<-lm(formula=formula,data=dat)
  testthat::expect_equal(res$Rsq[ind],summary(mod)$r.squared)
  testthat::expect_equal(res$AIC[ind],AIC(mod))
  testthat::expect_equal(res$BIC[ind],BIC(mod))
  
  ind<-8
  formula<-makeformula(lagAR=res$lagAR[ind],lag1=res$lag1[ind],lag2=res$lag2[ind])
  mod<-lm(formula=formula,data=dat)
  testthat::expect_equal(res$Rsq[ind],summary(mod)$r.squared)
  testthat::expect_equal(res$AIC[ind],AIC(mod))
  testthat::expect_equal(res$BIC[ind],BIC(mod))

  ind<-12
  formula<-makeformula(lagAR=res$lagAR[ind],lag1=res$lag1[ind],lag2=res$lag2[ind])
  mod<-lm(formula=formula,data=dat)
  testthat::expect_equal(res$Rsq[ind],summary(mod)$r.squared)
  testthat::expect_equal(res$AIC[ind],AIC(mod))
  testthat::expect_equal(res$BIC[ind],BIC(mod))
  
  #***test the function picklags_mat
  
  #make up some data
  set.seed(101)
  T<-100
  N<-3
  y<-matrix(rnorm(T*N),N,T)
  x1<-matrix(rnorm(T*N),N,T)
  x2<-matrix(rnorm(T*N),N,T)
  maxlags<-c(3,2,1)
  
  res<-picklags_mat(y,x1,x2,maxlags)  
  
  res1<-picklags_vect(y[1,],x1[1,],x2[1,],maxlags)
  res2<-picklags_vect(y[2,],x1[2,],x2[2,],maxlags)
  res3<-picklags_vect(y[3,],x1[3,],x2[3,],maxlags)
  testthat::expect_equal(res[,1:3],res1[,1:3])
  testthat::expect_equal(res[,1:3],res2[,1:3])
  testthat::expect_equal(res[,1:3],res3[,1:3])
  testthat::expect_equal(res[,4:6],(res1[,4:6]+res2[,4:6]+res3[,4:6])/3)
})
