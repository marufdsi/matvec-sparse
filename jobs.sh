qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v thread=$2,affinity=balance,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh

qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v thread=$2,affinity=scatter,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh

#qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v thread=$2,affinity=compact,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh

#qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
