
#qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$6:skylake -l mem=100GB -v thread=$1,row=$2,nnz=$3,run=$4,type=$5 -d $(pwd) balance.sh

for i in 1 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 225 250 275 300 400 500 600 700 800
do
	qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$5:skylake -l mem=100GB -v thread=$1,row=$((i * 1000000)),nnz=$((i * 1000000 * $2)),run=$3,type=$4 -d $(pwd) balance.sh
done

for i in 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 110 120 130 140 150 160 170 180 190 200 225 250 275 300 400 500 600 700 800
do
	qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$5:skylake -l mem=100GB -v thread=$1,row=$((i * 1000000)),nnz=$((i * 1000000 * $2)),run=$3,type=$4 -d $(pwd) balance.sh
done

#qsub -q copperhead -N "omp_spmv" -l nodes=1:ppn=$1:skylake -l mem=100GB -v thread=$1,row=$2,nnz=$3,run=$4,type=$5 -d $(pwd) balance.sh
