for size in 524288 1048576 2097152 4194304 
do
	for nz in 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288  
	do
		qsub -q copperhead -N "spmv_simple" -l nodes=$1:ppn=$2 -v thread=$3,run=$4,size=$size,nz=$nz -d $(pwd) runspMVOnServer.sh
	done
done
