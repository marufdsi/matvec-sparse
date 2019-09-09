qsub -q copperhead -N "seq_simple" -l nodes=1:ppn=1 -v parts=$2,matrix=$1 -d $(pwd) runSeqOnServer.sh

