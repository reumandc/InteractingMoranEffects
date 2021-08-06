#The main script, runs all other scripts, pulling necessary data from "Data" and putting all results
#in "Results"
#
#Before running this script, set the working directory of R to this directory.

rm(list=ls())
graphics.off()

#Clean the kelp data
source("Kelp_AnalysisStep_CleanData_SuperBasic.R")
source("Kelp_AnalysisStep_CleanData_MoreInvolved.R")

#Plot kelp time series for visual inspection
source("Kelp_AnalysisStep_PlotTimeSeries.R")

#Do linear models with kelp and related data to help demo the need for specific 
#lags - this also performs and saves the final level of data cleaning
source("Kelp_AnalysisStep_LinearModels.R")

#Do the spectra decomposition analysis for kelp
source("Kelp_AnalysisStep_SpectralDecomp.R")

#Clean the plankton data
