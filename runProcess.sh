qsub -q copperhead -N "mpi_simple" -l nodes=$1:ppn=$2 -v thread=$3,parts=$3,graph=$4 -d $(pwd) runOnServer.sh

