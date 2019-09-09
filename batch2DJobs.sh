#qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2:skylake -v procs=$3,input=$4 -d $(pwd) run2DSpMV.sh

qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v procs=$3,input=$4 -d $(pwd) run2DSpMV.sh
