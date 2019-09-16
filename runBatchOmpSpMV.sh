

for i in 36 32 18 16 9 8 5 4 3 2 1
do
	#qsub -q copperhead -N "omp_aff" -l nodes=1:ppn=$1:skylake -v file=$3,run=$2,affinity=$4,thread=$i,reserve=$1 -d $(pwd) runOmpSpMV.sh
	qsub -q copperhead -N "omp_aff" -l nodes=1:ppn=$1:skylake -v file=$3,run=$2,affinity=balanced,thread=$i,reserve=$1 -d $(pwd) runOmpSpMV.sh
#	qsub -q copperhead -N "omp_aff" -l nodes=1:ppn=$1:skylake -v file=$3,run=$2,affinity=scatter,thread=$i,reserve=$1 -d $(pwd) runOmpSpMV.sh
#	qsub -q copperhead -N "omp_aff" -l nodes=1:ppn=$1:skylake -v file=$3,run=$2,affinity=compact,thread=$i,reserve=$1 -d $(pwd) runOmpSpMV.sh
done

