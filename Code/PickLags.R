#This function picks lags in a linear regression setup relevant to the project. Read the specs and code for details.
#Just picks lags for each location, one by one, using picklags_vect, below, and then aggregates across locations. 
#
#Args
#y          An N by T matrix representing time series, N=number of locations, T=time series length. The population variable.
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
#Rsq        Rsq of the corresponding regression model, averaged across locations
#AIC        AIC of the corresponding regression model, averaged across locations
#BIC        BIC of the corresponding regression model, averaged across locations
#
#Notes
#When a given lag is included in the linear regression model, all smaller lags are automatically included.
#
picklags_mat<-function(y,x1,x2,maxlags)
{
  #compute results for each location and compute average results
  res_tot<-picklags_vect(y[1,],x1[1,],x2[1,],maxlags)
  for (counter in 2:(dim(y)[1]))
  {
    res_loc<-picklags_vect(y[counter,],x1[counter,],x2[counter,],maxlags)
    res_tot[,4:6]<-res_tot[,4:6]+res_loc[,4:6]
  }
  res<-res_tot
  res[,4:6]<-res[,4:6]/(dim(y)[1])
  
  #put in the delta AIC and BIC columns
  res$DeltaAIC<-res$AIC-min(res$AIC)
  res$DeltaBIC<-res$BIC-min(res$BIC)
  
  return(res)
}

#This function picks lags in a linear regression setup relevant to the project. Read the specs and code for details.
#This is for one location only, the above function iterates across locations.
#
#Args
#y          A vector of length T representing a time series, T=time series length. The population variable.
#x1         A vector of length T representing a time series. An environmental driver.
#x2         Another such
#maxlags    A vector of length 3. maxlags[1]>=0 is the max AR lag for the population variable, 0 meaning populations 
#             do not depend on past populations. maxlags[2]>=0 is the max lag for x1, 0 means only the current value
#             of x1 influences the populations. maxlags[3]>=0 is similar, but for x2.
#
#Output - a data frame with the following columns
#lagAR      The AR lag considered for the populations, >=0, meaning of 0 described above
#lag1       The lag for x1 considered, >=0, the meaning of 0 described above
#lag2       Similar but for x2
#Rsq        Rsq of the corresponding regression model
#AIC        AIC of the corresponding regression model
#BIC        BIC of the corresponding regression model
#
#Notes
#When a given lag is included in the linear regression model, all smaller lags are automatically included.
#
picklags_vect<-function(y,x1,x2,maxlags)
{
  #***modest error catching
  if ((length(y)!=length(x1)) || (length(y)!=length(x2)))
  {
    stop("Error in picklags_vect: y, x1 and x2 have to be vectors of the same length")
  }
  if (maxlags[1]<0 || maxlags[2]<0 || maxlags[2]<0)
  {
    stop("Error in picklags_vect: bad value for maxlags")
  }
  
  #***make the data frame for the regressions
  dat<-makedf(y,x1,x2,maxlags)

  #***iterate through all the regressions you want to do and do each and store results
  maxmaxlags<-max(maxlags)
  
  #make a storage receptacle for the results
  dim1res<-(maxlags[1]+1)*(maxlags[2]+1)*(maxlags[3]+1)
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
    for (lag1 in 0:maxlags[2])
    {
      for (lag2 in 0:maxlags[3])
      {
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
#the above function, separated out for unit testing purposes.
#
#Args
#y          A length-T time series, saved as a vector indexed by time, t. The population variable.
#x1         A length-T time series, saved as a vector indexed by time, t. An environmental driver.
#x2         Another such. 
#maxlags    A vector of length 3. maxlags[1]>=0 is the max AR lag for the population variable, 0 meaning populations 
#             do not depend on past populations. maxlags[2]>=-1 is the max lag for x1, -1 means x1 does not
#             influence the populations. maxlags[3]>=-1 is similar, but for x2.
#
#Output - The data frame necessary to do all the regressions performed by the function picklags_vect above.
#
#Notes - No error catching, not intended to be a user-facing function.
#
makedf<-function(y,x1,x2,maxlags)
{
  maxmaxlags<-max(maxlags)
  
  #put in the response variable
  dat<-data.frame(y=y[(1+maxmaxlags):length(y)])
  
  #put in the lagged versions of y
  if (maxlags[1]>0)
  {
    for (counter in 1:maxlags[1])
    {
      dat[,counter+1]<-as.vector(y[(1+maxmaxlags-counter):(length(y)-counter)])
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
      dat[,newcol]<-x1[(1+maxmaxlags-counter):(length(x1)-counter)]
      names(dat)[newcol]<-paste0("x1_l",counter)
      newcol<-newcol+1
    }
  }
  
  #put in the lagged versions of x2
  if (maxlags[3]>-1)
  {
    for (counter in 0:maxlags[3])
    {
      dat[,newcol]<-x2[(1+maxmaxlags-counter):(length(y)-counter)]
      names(dat)[newcol]<-paste0("x2_l",counter)
      newcol<-newcol+1
    }
  }
  
  return(dat)
}

#A function which constructs the regression formulas needs by the functions picklags_* above. 
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