#This function picks lags in a linear regression setup relevant to the project. Read the specs and code for details.
#
#Args
#y          An N by T matrix of time series, N=number of locations, T=time series length. The population variable.
#x1         An N by T matrix of time series. An environmental driver.
#x2         Another such
#maxlags    A vector of length 3. maxlags[1]>=0 is the max AR lag for the population variable, 0 meaning populations 
#             do not depend on past populations. maxlags[2]>=-1 is the max lag for x1, -1 means x1 does not
#             influence the populations. maxlags[3]>=-1 is similar, but for x2.
#
#Output - a data frame with the following columns
#lagAR      The AR lag considered for the populations, >=0, meaning of 0 described above
#lag1       The lag for x1 considered, >=-1, the meaning of -1 described above
#lag2       Similar but for x2
#Rsq        Rsq of the corresponding regression model
#AIC        AIC of the corresponding regression model
#BIC        BIC of the corresponding regression model
#
#Notes
#When a given lag is included in the linear regression model, all smaller lags are automatically included.
#
picklags<-function(y,x1,x2,maxlags)
{
  #***modest error catching
  if (any(dim(y)!=dim(x1)) || any(dim(y)!=dim(x2)))
  {
    stop("Error in picklags: y, x1 and x2 have to be the same dimensions")
  }
  if (maxlags[1]<0 || maxlags[2]<(-1) || maxlags[2]<(-1))
  {
    stop("Error in picklags: bad value for maxlags")
  }
  
  #***for convenience
  y<-t(y)
  x1<-t(x1)
  x2<-t(x2)
  
  #***make the data frame for the regressions
  dat<-makedf(y,x1,x2,maxlags)

  #***iterate through all the regressions you want to do and do each and store results
  maxmaxlags<-max(maxlags)
  
  #make a storage receptacle for the results
  dim1res<-(maxlags[1]+1)*(maxlags[2]+2)*(maxlags[3]+2)-1
  res<-data.frame(lagAR=NA*numeric(dim1res),
                  lag1=NA*numeric(dim1res),
                  lag2=NA*numeric(dim1res),
                  Rsq=NA*numeric(dim1res),
                  AIC=NA*numeric(dim1res),
                  BIC=NA*numeric(dim1res))
    
  #do the iteration
  currow<-1
  for (lagAR in 0:maxlags[1])
  {
    for (lag1 in -1:maxlags[2])
    {
      for (lag2 in -1:maxlags[3])
      {
        if(lagAR==0 && lag1==-1 && lag2==-1) #the case with no predictors not considered
        {
          next;
        }
        
        #make the formula for the call to lm
        formula<-makeformula(lagAR,lag1,lag2)
          
        #call lm
        mod<-lm(formula=formula,data=dat)
          
        #save results in the appropriate place in res
        res[currow,1]<-lagAR
        res[currow,2]<-lag1
        res[currow,3]<-lag2
        res[currow,4]<-summary(mod)$r.squared
        res[currow,5]<-AIC(mod)
        res[currow,6]<-BIC(mod)
        
        currow<-currow+1
      }
    }
  }
  
  #return results
  return(res)
}

#Makes the data frame for the regressions in the above function. This is a utility function serving
#the above functiuon, separated out for unit testing purposes.
#
#Args
#y          A T by N matrix of time series, N=number of locations, T=time series length. The population variable.
#x1         A T by N matrix of time series. An environmental driver.
#x2         Another such. Note these are transposed from what is expected in the function picklags above.
#maxlags    A vector of length 3. maxlags[1]>=0 is the max AR lag for the population variable, 0 meaning populations 
#             do not depend on past populations. maxlags[2]>=-1 is the max lag for x1, -1 means x1 does not
#             influence the populations. maxlags[3]>=-1 is similar, but for x2.
#
#Output - The data frame necessary to do all the regressions performed by the function picklags above.
#
#Notes - No error catching, not intended to be a user-facing function.
#
makedf<-function(y,x1,x2,maxlags)
{
  maxmaxlags<-max(maxlags)
  
  #put in the response variable
  dat<-data.frame(y=as.vector(y[1:(dim(y)[1]-maxmaxlags),]))
  
  #put in the lagged versions of y
  if (maxlags[1]>0)
  {
    for (counter in 1:maxlags[1])
    {
      dat[,counter+1]<-as.vector(y[(1+counter):(dim(y)[1]-maxmaxlags+counter),])
      names(dat)[counter+1]<-paste0("y_l",counter)
    }
  }
  
  #put in the lagged versions of x1
  newcol<-1+maxlags[1]+1 #the response variable, plus the maxlags[1] AR lags, plus 1 for the 
                         #position of where the new column will be 
  if (maxlags[2]>-1)
  {
    for (counter in 0:maxlags[2])
    {
      dat[,newcol]<-as.vector(x1[(1+counter):(dim(x1)[1]-maxmaxlags+counter),])
      names(dat)[newcol]<-paste0("x1_l",counter)
      newcol<-newcol+1
    }
  }
  
  #put in the lagged versions of x2
  if (maxlags[3]>-1)
  {
    for (counter in 0:maxlags[3])
    {
      dat[,newcol]<-as.vector(x2[(1+counter):(dim(x2)[1]-maxmaxlags+counter),])
      names(dat)[newcol]<-paste0("x2_l",counter)
      newcol<-newcol+1
    }
  }
  
  return(dat)
}

#A function which constructs the regression formulas needs by the function picklags above. 
#Separated out for unit testing purposes.
#
#Args
#lagAR      The autoregressive lags used are 1 up to this value, if it is 0 no AR lags are used
#lag1       The x1 lags used are 0 up to this value, if it is -1, no x1 lags are used
#lag2       The x2 lags used are 0 up to this value, if it is -1, no x2 lags are used
#
#Output - an object of class formula for the regression
#
#Notes
#No intercept in the regression
#No error catching, not a user-facing function
#Assumed lagAR>=0 an integer
#Assumed lag1>=-1 an integer
#Assumed lag2>=-1 an integer
#Assumed we are not in the case lagAR==0 and lag1==-1 and lag2==-1
#
makeformula<-function(lagAR,lag1,lag2)
{
  #the response variable
  formula<-"y~"
  
  #the AR lags
  if (lagAR>0)
  {
    for (counter in 1:lagAR)
    {
      if (counter==1)
      {
        formula<-paste0(formula,"y_l1")
      }else
      {
        formula<-paste0(formula,"+y_l",counter)
      }
    }
  }
  
  #the x1 lags
  if (lag1>-1)
  {
    for (counter in 0:lag1)
    {
      if (lagAR==0 && counter==0)
      {
        formula<-paste0(formula,"x1_l0") 
      }else
      {
        formula<-paste0(formula,"+x1_l",counter)
      }
    }
  }
  
  #the x2 lags
  if (lag2>-1)
  {
    for (counter in 0:lag2)
    {
      if (lagAR==0 && lag1==-1 && counter==0)
      {
        formula<-paste0(formula,"x2_l0") 
      }else
      {
        formula<-paste0(formula,"+x2_l",counter)
      }
    }
  }
  
  #no intercept term
  formula<-paste0(formula,"-1")
  
  #convert to formula class and return
  formula<-as.formula(formula)
  return(formula)
}