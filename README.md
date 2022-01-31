# How environmental drivers of synchrony interact - Introduction to the repository of analyses supporting the paper

Daniel C. Reuman, University of Kansas, reuman@ku.edu

Max Castorani, University of Virginia

Kyle Cavanaugh, UCLA

Lawrence W. Sheppard, Marine Biological Association of the United Kingdom

Jonathan Walter, University of Virginia

Tom Bell, UCSB and Woods Hole Oceanographic Institution

## Introduction

This repository can be used to reproduce the analyses behind the paper "How environmental drivers of synchrony interact" and to recompile the latex of the paper. Data are also included in the repository. 

## How to compile

If you want to run the unit tests of the function used, make your R working directory equal to the Code directory of the repository and run `testthat::test_dir(".")`. This should take less than 1 minute.

To run the codes that reproduce the analyses of the paper, make your R working directory equal to the `Code` directory of the repository and run `source("main.R")`. If all dependencies are in place, this will pull data from the `Data` directory and create results and then put them in the `Results` directory. All results are regenerated except Fig. 2, which was produced separately by Max Castorani on his machine, based on data exported by the script `Kelp_AnalysisStep_MakeMapFigDat.R` which is invoked by `main.R`. It was not considered necessary to embed Max's scripts into this workflow because Fig. 2 is purely descriptive, intended solely to introduce the kelp system; main results are not based on Fig. 2. The run will take some time the first time around because of use of the `checkpoint` package (see below). Subsequent runs (if any) should be faster because `checkpoint` will already have installed all the necessary R packages.

To recompile the main text and supporting information of the paper, knit the Latex/Sweave files `Paper.Rnw` and `SupMat.Rnw`, which are in the `Paper` directory. These files will import results from the `Results` directory and will create two pdfs in the `Paper` directory.

## Dependencies

### Overview of dependencies

R, R studio, checkpoint package, a setup that can compile Sweave documents. 

### Dependencies on R, R studio

I used R version 3.6.3 running on Ubuntu 16.04, and R studio version 1.1.423.

### Dependencies on the R `checkpoint` package

This codes uses the R `checkpoint` package. This is set up in the master file `main.R`, which contains a line of code specifying a date.

`checkpoint::checkpoint(snapshot_date="2022-01-01",r_version=getRversion(),checkpoint_location=getwd(), scan_now = TRUE)` 

The `checkpoint` package then automatically scans through other files looking for other required R packages. It then downloads and installs the newest versions of those packages available on the given date. This helps ensure that re-compiling the document uses exactly the same code that was originally used. This can take some time on first run (as mentioned above) but it is faster on subsequent runs because the packages are already installed. This also means that R package dependencies should only be the `checkpoint` package, since that package should scan for other packages and install them locally. Quite a few MB disk space are used (80-100). The version `checkpoint` that I used was 1.0.2. 

### Dependencies on a setup for compiling Seave documents

Set your R studio up so there is a button near the top that compiles .Rnw files into pdfs.

### Additional dependencies

The compilation process was tested by Reuman on Ubuntu 16.04 using R version 3.6.3 and R studio version 1.1.423. It has not been tested on other machines. We have endeavored to list all dependencies we can think of above, but we have only compiled on Reuman's machine, so we cannot guarantee that additional dependencies will not also be needed on other machines. This repository is intended to record a workflow, and is not designed or tested for distribution and wide use on multiple machines. It is not guaranteed to work on the first try without any hand-holding on arbitrary computing setups.

## Intermediate files:

Knitting the Sweave code automatically produces a lot of 'intermediate' files. Files ending in `.tex` are the converted documents from `.Rnw` including all the R code output and the rest (files ending `.log`, `.aux`, `.lof`, `.lot`, `.toc`  and `.out` ) are intermediate files that `pdflatex` uses to keep track of various parts of the document. Some of these can be useful for diagnosing problems, if any. 

## Acknowlegements

This study was partly funded by the U.S. National Science Foundation (NSF) through linked NSF-OCE awards 
2023555, 2023523, 2140335, and 2023474 to M.C.N.C, K.C.C., T.W.B, and D.C.R., respectively; and by NSF-Math Bio award 
1714195 to D.C.R.; and by support to D.C.R. from the James S. McDonnell Foundation and the California Department of Fish
and Wildlife Delta Science Program. This project used data developed through the Santa Barbara Coastal 
Long Term Ecological Research project, funded through NSF-OCE award 1831937. The authors thank Vadim Karatayev, 
Maowei Liang, Kyle Emery, Nat Coombs, Adeola Adeboje and Ethan Kadiyala for helpful discussions.
%***DAN: If you end up adding any of these people to the author list, take them off here.
Any opinions, findings, and conclusions or recommendations expressed in this material are those of 
the authors and do not necessarily reflect the views of the National Science Foundation or the other funders. 












