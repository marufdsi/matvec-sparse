


#qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v procs=$3 -d $(pwd) run1DSpMV.sh


qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=nlpkkt160,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=nlpkkt240,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=cage15,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=G3_circuit,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=kkt_power,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=ecology1,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=uk-2002,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=asia,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=NLR,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=coPapersDBLP,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=coPapersCiteseer,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=copter2,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=333SP,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=af_shell10,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=AS365,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=europe,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=M6,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=nlpkkt200,procs=$3 -d $(pwd) run1DSpMV.sh
qsub -q copperhead -N "omp_aff" -l nodes=$1:ppn=$2 -v input=thermal2,procs=$3 -d $(pwd) run1DSpMV.sh
