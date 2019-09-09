for size in 524288 1048576 2097152
do
	col = size*$3
  	for nz in 131072 262144
        do
          	qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$size,col=$col,nz=$nz -d $(pwd) runSpMVModel.sh
        done
done

for size in 4194304
do
	col = size*$3
  	for nz in 524288 1048576 2097152 4194304
        do
          	qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$size,col=$col,nz=$nz -d $(pwd) runSpMVModel.sh
        done
done

for size in 134217728
do
	col = size*$3
  	for nz in 8388608 16777216 33554432 67108864 67108864
        do
          	qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=$size,col=$col,nz=$nz -d $(pwd) runSpMVModel.sh
        done
done

qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=32768,col=$(expr 32768 '*' $3),nz=65536 -d $(pwd) runSpMVModel.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=65536,col=$(expr 65536 '*' $3),nz=131072 -d $(pwd) runSpMVModel.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=131072,col=$(expr 131072 '*' $3),nz=262144 -d $(pwd) runSpMVModel.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=262144,col=$(expr 262144 '*' $3),nz=524288 -d $(pwd) runSpMVModel.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=524288,col=$(expr 524288 '*' $3),nz=1048576 -d $(pwd) runSpMVModel.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=1048576,col=$(expr 1048576 '*' $3),nz=2097152 -d $(pwd) runSpMVModel.sh
qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=2097152,col=$(expr 2097152 '*' $3),nz=4194304 -d $(pwd) runSpMVModel.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=4194304,col=$(expr 4194304 '*' $3),nz=8388608 -d $(pwd) runSpMVModel.sh
#qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,row=8388608,col=$(expr 8388608 '*' $3),nz=16777216 -d $(pwd) runSpMVModel.sh
