#!/bin/sh
#PBS -l ncpus=32,mpiprocs=32,mem=500gb,host=hpc_c11
#PBS -M fujimori.shinichiro.8a@kyoto-u.ac.jp
#PBS -m abe
#PBS -N AIMHub-vis
#PBS -e ../../../../output/jobreport
#PBS -o ../../../../output/jobreport
export OMP_NUM_THREADS=32      # 環境変数の設定
. ~/.bashrc
echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
echo `pwd`
timestamp=`date +"%m-%d-%y-%H-%M-%S"`
echo $timestamp > ../../../../output/jobreport/psublog_${PBS_JOBID}.txt 
cat ./settings/vis.sh >> ../../../../output/jobreport/psublog_${PBS_JOBID}.txt
bash AR6.sh >> ../../../../output/jobreport/psublog_${PBS_JOBID}.txt
