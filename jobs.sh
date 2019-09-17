#qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v thread=$2,affinity=balanced,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh

#qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v thread=$2,affinity=scatter,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh

#qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v thread=$2,affinity=compact,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh

#qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh


qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/nlpkkt160_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/nlpkkt240_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/cage15_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/G3_circuit_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/kkt_power_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/ecology1_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/uk-2002_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/asia_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/NLR_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/coPapersDBLP_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/coPapersCiteseer_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/copter2_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/333SP_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/af_shell10_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/AS365_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/europe_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/M6_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/nlpkkt200_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh
qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -v file=input/thermal2_original.mtx,thread=$2,affinity=$4,reserve=$1,run=$3 -d $(pwd) runGraphOnOmpSpMV.sh

