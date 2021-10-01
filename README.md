# AIMCGEVisualization
Visualizing AIMCGE output

Date: 2018/7/17 updated in 2020/6/1 (2021/10/01)

Name: Shinichiro Fujimori

- Overview

This tool is for the standarized visualization of AIMCGE(AIMHub) output.


- Input

In /data/ folder, there are input files


- region.txt : region set
- varlist.txt : variable list (using AIM IIASA coding system)
- scenariomap.map :scenario list and its names
- specify whether AIM/Enduse results are included in the analysis or not by "enduseflag"

In /modeloutput/ folder, the model output should be included.
- global_17_emf.gdx : input file

- Location
This tool can be used in the same folder as AIM/CGE,
but you can also use it in the independent folder which requires folder specification to copy the above CGE results GDX files.

- Output
outout of this tool is in /output directory where you see regional names and "merge". The latter combine the regional resutls and the coverage of the regions should be specified in data/region.txt
 
