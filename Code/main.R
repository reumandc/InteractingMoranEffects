#The main script, runs all other scripts, pulling necessary data from "Data" and putting all results
#in "Results"
#
#Before running this script, set the working directory of R to this directory.

rm(list=ls())
graphics.off()

#Clean the data
source("AnalysisStep_CleanData_SuperBasic.R")
source("AnalysisStep_CleanData_MoreInvolved.R")

#plot time series for visual inspection
source("AnalysisStep_PlotTimeSeries.R")

#Do linear models to help demo the need for specific lags - this also performs and saves the final level of data cleaning
source("AnalysisStep_LinearModels.R")

#Do the spectra decomposition analysis
source("AnalysisStep_SpectralDecomp.R")

