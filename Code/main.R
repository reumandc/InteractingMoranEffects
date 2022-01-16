#The main script, runs all other scripts, pulling necessary data from "Data" and putting all results
#in "Results"
#
#Before running this script, set the working directory of R to this directory.

rm(list=ls())
graphics.off()

#Develop the theoretical case studies
source("Theory_CaseStudyAB.R")
source("Theory_CaseStudyC.R")
source("Theory_IntroExample.R")
source("Theory_PedagogFig.R")

#Clean the kelp data
source("Kelp_AnalysisStep_CleanData_SuperBasic.R")
source("Kelp_AnalysisStep_CleanData_MoreInvolved.R")

#Plot kelp time series for visual inspection
#source("Kelp_AnalysisStep_PlotTimeSeries.R")

#Final cleaning
source("Kelp_AnalysisStep_CleanData_Final.R")

#Make the map fig, to help orient the reader of the paper to the dataset
source("Kelp_AnalysisStep_MakeMapFig.R")

#Do the spectra decomposition analysis for kelp
source("Kelp_AnalysisStep_SpectralDecomp.R")


