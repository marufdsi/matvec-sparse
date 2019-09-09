

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
do
	#qsub -q copperhead -N "omp_aff" -l nodes=1:ppn=$1:skylake -v file=$3,run=$2,affinity=$4,thread=$i,reserve=$1 -d $(pwd) runOmpSpMV.sh
	qsub -q copperhead -N "omp_aff" -l nodes=1:ppn=$1:skylake -v file=$3,run=$2,affinity=balance,thread=$i,reserve=$1 -d $(pwd) runOmpSpMV.sh
	qsub -q copperhead -N "omp_aff" -l nodes=1:ppn=$1:skylake -v file=$3,run=$2,affinity=scatter,thread=$i,reserve=$1 -d $(pwd) runOmpSpMV.sh
	qsub -q copperhead -N "omp_aff" -l nodes=1:ppn=$1:skylake -v file=$3,run=$2,affinity=compact,thread=$i,reserve=$1 -d $(pwd) runOmpSpMV.sh
done

