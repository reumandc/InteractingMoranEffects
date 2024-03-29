# How environmental drivers of synchrony interact - Introduction to the repository of analyses supporting the paper

Daniel C. Reuman, University of Kansas, reuman@ku.edu  
Max Castorani, University of Virginia  
Kyle Cavanaugh, UCLA  
Lawrence W. Sheppard, Marine Biological Association of the United Kingdom  
Jonathan Walter, University of Virginia  
Tom Bell, UCSB and Woods Hole Oceanographic Institution  

## Introduction

This repository can be used to reproduce the analyses behind the paper "How environmental drivers of synchrony interact" and to recompile the latex of the paper. Data are also included in the repository. As of June 2023, the paper has now been accepted for publication
in the journal Ecography. 

## How to compile

If you want to run the unit tests of the function used, make your R working directory equal to the Code directory of the repository and run `testthat::test_dir(".")`. This should take less than 1 minute. You may have to install an R version or some R packages
that you don't have to get these tests to pass. Keep in mind these tests were run when the paper was written, with versions 
R and R studio and R packages which were current around then (see below), so no guarantees if a lot of time has passed since then.
You can always run the code on the (now) old version of R and the packages and it should work then, but you'll have to set
up a (now) old version of R with (now) old packages yourself. Yes, I have heard of the checkpoint package, Microsoft has now
abandoned MRAN so this is not a way toward reproducibility anymore.

To run the codes that reproduce the analyses of the paper, make your R working directory equal to the `Code` directory of 
the repository 
and run `source("main.R")`. If all dependencies are in place (see below), this will pull data from the `Data` directory and create 
results and then put them in the `Results` directory. All results are regenerated except Fig. 3, which was produced separately by Max 
Castorani on his machine, based on data exported by the script `Kelp_AnalysisStep_MakeMapFigDat.R` which is invoked by `main.R`. It was 
not considered necessary to embed Max's scripts into this workflow because Fig. 3 is purely descriptive, intended solely to introduce 
the 
kelp system; main results are not based on Fig. 3. 
 
To recompile the main text and supporting information of the paper, knit the Latex/Sweave files `Paper.Rnw` and `SupMat.Rnw`, which are 
in the `Paper` directory. These files will import results from the `Results` directory and will create two pdfs in the `Paper` 
directory. 
You may have to do several knits, alternating between the two files, to get all the links to work. 

## Dependencies

### Overview of dependencies

R, R studio, a setup that can compile Sweave documents. The right versions of everything may be required - see below.

### Dependencies on R, R studio

I used R version 4.3.0 running on Ubuntu 18.04, and R studio version 2023.03.0+386. 

### Dependencies on a setup for compiling Seave documents

Set your R studio up so there is a button near the top that compiles .Rnw files into pdfs.

### Package versions

The packages I used were current around the time of R version 4.3.0. Specifically:
matrixcalc, version 1.0-6;
wsyn, version 1.0.4;
ncf, version 1.3-2;
latex2exp, version 0.9.6; 
testthat, version 3.1.9;
mvtnorm, version 1.2-2;
shape, version 1.4.6;
graphics, version 4.3.0;
stats, version 4.3.0;
parallel, version 4.3.0.

You may have to install the versions of R, R studio, and the packages listed above to get the code to
run if a lot of time has passed since publication. 

### Additional dependencies

The code and paper compilation process was tested by Reuman on Ubuntu 18.04 using R version 4.3.0 and R studio version 2023.03.0+386. 
It has not been tested on other machines. We have endeavored to list all dependencies we can think of above, but we have 
only run the code on Reuman's machine, so we cannot guarantee that additional dependencies will not also be needed on other 
machines. This repository is intended to record a workflow, and is not designed or tested for distribution and wide use on 
multiple machines. It is not guaranteed to work on the first try without any hand-holding on arbitrary computing setups.

## Intermediate files:

Knitting the Sweave code automatically produces a lot of 'intermediate' files. Files ending in `.tex` are the converted documents 
from `.Rnw` including all the R code output and the rest (files ending `.log`, `.aux`, `.lof`, `.lot`, `.toc`  and `.out` ) 
are intermediate files that `pdflatex` uses to keep track of various parts of the document. Some of these can be useful for 
diagnosing problems, if any. 

## Acknowlegements

This study was partly funded by the U.S. National Science Foundation (NSF) through linked NSF-OCE awards 
2023555, 2023523, 2140335, and 2023474 to M.C.N.C, K.C.C., T.W.B, and D.C.R., respectively; and by NSF-Math Bio award 
1714195 to D.C.R.; and by support to D.C.R. from the James S. McDonnell Foundation, the Humboldt Foundation, 
and the California Department of Fish and Wildlife Delta Science Program. This project used data developed 
through the Santa Barbara Coastal Long Term Ecological Research project, funded through NSF-OCE award 1831937. 
The authors thank Vadim Karatayev, Maowei Liang, Kyle Emery, Nat Coombs, Adeola Adeboje and Ethan Kadiyala for 
helpful discussions. Any opinions, findings, and conclusions or recommendations expressed in this material are those of 
the authors and do not necessarily reflect the views of the National Science Foundation or the other funders. 












