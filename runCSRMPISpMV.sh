#! /bin/bash

qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -l mem=300GB -v thread=$3,run=$4,row=$5,col=$(expr $3 '*' $5),nz=$6 -d $(pwd) runCSRSpMV.sh

#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=1024 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=512 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=256 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=128 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=64 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=32 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=16 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=8 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=2048 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=4096 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=8192 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=16384 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=32768 -d $(pwd) runCSRSpMV.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$5,nz=65536 -d $(pwd) runCSRSpMV.sh

