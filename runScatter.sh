qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$6:skylake -l walltime=10:00:00 -v thread=$1,row=$2,nnz=$3,run=$4,type=$5 -d $(pwd) scatter.sh
#qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -l mem=100GB -v thread=$1,row=$2,nnz=$3,run=$4,type=$5 -d $(pwd) scatter.sh

