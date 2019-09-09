qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v procs=$3 -d $(pwd) run1DSpMV.sh
