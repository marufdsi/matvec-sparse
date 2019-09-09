

qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$6:skylake -l mem=100GB -v thread=$1,row=$2,nnz=$3,run=$4,type=$5 -d $(pwd) compact.sh



#qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -l mem=100GB -v thread=$1,row=$2,nnz=$3,run=$4,type=$5 -d $(pwd) compact.sh

