#!/bin/bash -x
#. ~/.bashrc
export gams_sys_dir=`which gams --skip-alias|xargs dirname`
NCPU=30
#arguments are: 1:gams sys,2:number of CPU, 3:visualizaton scenario name specification auto or not, 4:file location, 5:submodule switch, 
#               6:enduse iteration switch, 7: GDX file name, 8: region code/global, 9: AR6 consideration, 10: IntTool
echo   R --vanilla --slave --args ${gams_sys_dir} ${NCPU} < AR6.r
R --vanilla --slave --args ${gams_sys_dir} ${NCPU} < AR6.r
cd ../data
tar -zcvf ./AR6Slected.tar.gz ./AR6Slected.csv
cd ../prog