#!/usr/bin/bash
#SBATCH -p general # partition (queue)
#SBATCH -N 1                      # number of nodes
#SBATCH -J bowtie2
#SBATCH -n 4                      # number of cores
#SBATCH --mem=8G
#SBATCH --qos=general
#SBATCH -o bowtie2.%N.%A_%a.out        # STDOUT
#SBATCH --mail-type=END,FAIL      # notifications for job done & fail
#SBATCH --mail-user=cera.fisher@uconn.edu # send-to address

## Set up needful environment variables
source ~/.profile
LIBLIST=$1
module load bowtie2/2.2.9

## Read the libs file to get the array of library prefixes.
mapfile -t LIB < $LIBLIST

CURRENT_LIB=${LIB[$SLURM_ARRAY_TASK_ID]}
LIBEXT=`grep $CURRENT_LIB 2020_library_basenames.txt | cut -f 2`
RNADB=`grep $CURRENT_LIB 2020_library_information.txt | cut -f 6`

cp $LABS_DIR/AssembledReads/*$CURRENT_LIB*P1.trimmed2.fastq.gz input/$CURRENT_LIB\.$LIBEXT\.R1.fastq.gz
cp $LABS_DIR/AssembledReads/*$CURRENT_LIB*P2.trimmed2.fastq.gz input/$CURRENT_LIB\.$LIBEXT\.R2.fastq.gz
cp $LABS_DIR/AssembledReads/*$CURRENT_LIB*U1.trimmed2.fastq.gz input/$CURRENT_LIB\.$LIBEXT\.U1.fastq.gz
cp $LABS_DIR/AssembledReads/*$CURRENT_LIB*U2.trimmed2.fastq.gz input/$CURRENT_LIB\.$LIBEXT\.U2.fastq.gz

R1COUNT=`zgrep -c '^+$' input/$CURRENT_LIB\.$LIBEXT\.R1.fastq.gz`
U1COUNT=`zgrep -c '^+$' input/$CURRENT_LIB\.$LIBEXT\.U1.fastq.gz`
U2COUNT=`zgrep -c '^+$' input/$CURRENT_LIB\.$LIBEXT\.U2.fastq.gz`

echo "Input R1 Reads: $R1COUNT"
echo "Input U1 Reads: $U1COUNT"
echo "Input U2 Reads: $U2COUNT"

bash bowtie2_func.sh $CURRENT_LIB\.$LIBEXT $RNADB ERCC92
wait

#EJCF002_N01.EbinL5T_Leg.ERCC92.bam

if [ -s output/$CURRENT_LIB\.$LIBEXT\.ERCC92.bam ]; then
    rm -f output/$CURRENT_LIB\.*.bam
fi

mv bowtie2.*.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID\.out filter_logs/.

