#The main script, runs all other scripts, pulling necessary data from "Data" and putting all results
#in "Results"
#
#Before running this script, set the working directory of R to this directory.

#Clean the data
source("CleanData_SuperBasic.R")

#
source("Approach_ModellingSpectra_General.R")

