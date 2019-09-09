for size in 524288 1048576 2097152
do
	for nz in 131072 262144
	do
		qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=$size,nz=$nz -d $(pwd) runspMVOnServer.sh
	done
done

for size in 4194304                                                     
do
       	for nz in 524288 1048576 2097152 4194304
        do
                qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=$size,nz=$nz -d $(pwd) runspMVOnServer.sh
        done
done

for size in 134217728
do
       	for nz in 8388608 16777216 33554432
	do
		qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=$size,nz=$nz -d $(pwd) runspMVOnServer.sh
	done
done
