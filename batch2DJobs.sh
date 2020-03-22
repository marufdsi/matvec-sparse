#qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2:skylake -v procs=$3,input=$4 -d $(pwd) run2DSpMV.sh

#qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v procs=$3,input=$4 -d $(pwd) run2DSpMV.sh

qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v procs=$3 -d $(pwd) run2DSpMV.sh


#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=nlpkkt160,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=nlpkkt240,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=cage15,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=G3_circuit,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=kkt_power,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=ecology1,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=uk-2002,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=asia,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=NLR,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=coPapersDBLP,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=coPapersCiteseer,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=copter2,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=333SP,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=af_shell10,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=AS365,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=europe,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=M6,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=nlpkkt200,procs=$3 -d $(pwd) run2DSpMV.sh
#qsub -q copperhead -N "mpi_2D" -l nodes=$1:ppn=$2:skylake -v input=thermal2,procs=$3 -d $(pwd) run2DSpMV.sh

