#The main script, runs all other scripts, pulling necessary data from "Data" and putting all results
#in "Results"
#
#Before running this script, set the working directory of R to this directory.

rm(list=ls())
graphics.off()

#***
#checkpoint, to help ensure it continues to work into the future
#***

#print("Running checkpoint, this may take a while the first time...")
#if (!dir.exists("./.checkpoint/")){
#  dir.create("./.checkpoint/")
#}
#checkpoint::checkpoint(snapshot_date="2023-06-01",r_version=getRversion(),checkpoint_location=getwd(),
#                       scan_now = TRUE) 

#***
#Develop the theoretical case studies
#***

print("Theory case studies A and B...")
source("Theory_CaseStudyAB.R")

print("Theory case study C...")
source("Theory_CaseStudyC.R")

print("Theory example in intro...")
source("Theory_IntroExample.R")

print("Pedagogical figures...")
source("Theory_PedagogFig_Alt.R")

#***
#Clean the data
#***

print("Super basic data cleaning...")
source("Kelp_AnalysisStep_CleanData_SuperBasic.R")

print("More involved data cleaning...")
source("Kelp_AnalysisStep_CleanData_MoreInvolved.R")

#print("Plot kelp time series for visual inspection...")
#source("Kelp_AnalysisStep_PlotTimeSeries.R")

print("Final cleaning...")
source("Kelp_AnalysisStep_CleanData_Final.R")

#***
#Make the map figure, to help orient the reader of the paper to the dataset 
#***

print("Mapping the system...")
source("Kelp_AnalysisStep_MakeMapFigDat.R") #preps the data for Max Castorani code,
#run separately on another machine to produce Fig. 2, which is purely descriptive

#***
#Do the spectra decomposition analysis for kelp
#***

print("Spectral analysis...")
source("Kelp_AnalysisStep_SpectralDecomp.R")
print("Done!")

checkpoint::uncheckpoint()
