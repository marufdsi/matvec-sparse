#! /bin/bash

qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=32768,nz=65536 -d $(pwd) runspMVOnServer.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=65536,nz=131072 -d $(pwd) runspMVOnServer.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=131072,nz=262144 -d $(pwd) runspMVOnServer.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=262144,nz=524288 -d $(pwd) runspMVOnServer.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=524288,nz=1048576 -d $(pwd) runspMVOnServer.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=1048576,nz=2097152 -d $(pwd) runspMVOnServer.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=2097152,nz=4194304 -d $(pwd) runspMVOnServer.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=4194304,nz=8388608 -d $(pwd) runspMVOnServer.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=8388608,nz=16777216 -d $(pwd) runspMVOnServer.sh
