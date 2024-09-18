# AIMCGEVisualization
Visualizing AIMCGE output

Date: 2018/7/17 updated in 2020/6/1 (2021/10/01)

Name: Shinichiro Fujimori

# Overview

This tool is for the standarized visualization of AIMCGE(AIMHub) output.


# Input
## Specification by files

In /data/ folder, there are input files

- region.txt : region set
- varlist.txt : variable list (using AIM IIASA coding system)
- scenariomap.map :scenario list and its names
- scenariomap2.map :scenario list and its names for Endduse (loading process)
- EnduseScenarioMap.map :scenario list and its names for Endduse (final presetation)
- specify whether AIM/Enduse results are included in the analysis or not by "enduseflag"

In /modeloutput/ folder, the model output should be included. This is automated by the program but may be better double check.
- global_17_iamc.gdx : input file

# Spepecification by parameters

In the main program there is a location to specifies key parameters which tag is below  
  #---------------switches to specify the run condition -----  
filename <- "global_17" # filename should be "global_17","CHN","JPN"....  
enduseflag <- 0   # If you would like to display AIM/Enduse outputs, make this parameter 1 otherwise 0.  
dirCGEoutput <-"../../anls_output/iiasa_database/gdx/"  # directory where the CGE output is located   
parallelmode <- 1 #Switch for parallel process. if you would like to use multi-processors assign 1 otherwise 0.  
threadsnum <- min(floor(availableCores()/2),24)  
r2ppt <- 0 #Switch for ppt export. if you would like to export as ppt then assign 1 otherwise 0.  
mergecolnum <- 6 #merge figure facet number of columns  
CGEgdxcopy <- 0 # if you would like to copy and store the CGE IAMC template file make this parameter 1, otherwise 0.  


# Location:
This tool can be used in the same folder as AIM/CGE,
but you can also use it in the independent folder which requires folder specification to copy the above CGE results GDX files.

# Output:
outout of this tool is in /output directory where you see regional names and "merge". The latter combine the regional resutls and the coverage of the regions should be specified in data/region.txt
 

# Notes on changing the code
## Add new area figure

The following files should be edited

- ../data/Areafigureorder.txt
- ../data/Area.map

The first file identifies the list of indicators, units and y axis label.
The second file describes what indicators are assigned to each area plot, order of them and color schemes.


#Required R packages
- reshape2
- cowplot
- ggplot2
- RColorBrewer
- dplyr
- sp
- maptools
- maps
- ggradar
- fmsb
- tidyr
- stringr
- rJava
- Rcpp
- ReporteRsjars
- ReporteRs
- xlsx
- officer
- furrr
- purrr
- progressr
- gdxrrw


  liblist <- c("reshape2","cowplot","ggplot2","RColorBrewer","dplyr","sp","maps","fmsb","tidyr","stringr","rJava","Rcpp","ReporteRs","xlsx","officer","furrr","purrr","progressr","ncdf4","openxlsx","gdata","haven","reldist","foreign","plyr","tidyverse","readxl","raster","foreach","doParallel","gdxrrw")
for(j in liblist){
  eval(parse(text=paste0("library(",j,")")))
}

"maptools",ggradar
ReporteRsjars