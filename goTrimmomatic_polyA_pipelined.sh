#!/bin/bash
#
#SBATCH -p general # partition (queue)
#SBATCH -N 1                      # number of nodes
#SBATCH -J runTrim
#SBATCH -n 4                      # number of cores
#SBATCH --mem=8G
#SBATCH --qos=general
#SBATCH -o trimmomatic.%N.%A_%a.out        # STDOUT
#SBATCH -e trimmomatic.%N.%A_%a.err        # STDERR
#SBATCH --mail-type=END,FAIL      # notifications for job done & fail
#SBATCH --mail-user=cera.fisher@uconn.edu # send-to address

## Make these directories in your workdir before running:
# polyA_output/
# polya_logs/

mapfile -t LIB < libs.txt

CURRENT_LIB=${LIB[$SLURM_ARRAY_TASK_ID]}
LIBEXT=`grep "$CURRENT_LIB" 2020_library_basenames.txt | cut -f 2`

READDIR=output/

wheretrim=/isg/shared/apps/Trimmomatic/0.36
LEFT=$CURRENT_LIB.P1.trimmed.fastq.gz
RIGHT=$CURRENT_LIB.P1.trimmed.fastq.gz
#cp EJCF002_N11_S221_L003_R1_001.fastq.gz ~/bio/Membracids/ReadProcessing/input/.

if [ -s $READDIR/$LEFT ];
then
    module load Trimmomatic
    java -jar $wheretrim/trimmomatic-0.36.jar PE -threads 4 $READDIR/$LEFT $READDIR/$RIGHT polyA_output/$CURRENT_LIB_$LIBEXT.P1.trimmed2.fastq.gz polyA_output/$CURRENT_LIB.U1.trimmed2.fastq.gz polyA_output/$CURRENT_LIB_$LIBEXT.P2.trimmed2.fastq.gz polyA_output/$CURRENT_LIB.U2.trimmed2.fastq.gz ILLUMINACLIP:adapters/polyAtails.fa:2:30:10;
fi
### These will be trimmed2 because I can't think of a better less cumbersome name to mark what part
### of the process they've gone through

if [ -s polyA_output/$CURRENT_LIB_$LIBEXT.P1.trimmed2.fastq.gz ];
then
   touch $CURRENT_LIB.polyA.done;
   echo "$CURRENT_LIB   $CURRENT_LIB    polyA_output/$CURRENT_LIB_$LIBEXT.P1.trimmed2.fastq.gz  polyA_output/$CURRENT_LIB_$LIBEXT.P2.trimmed2.fastq.gz" > $CURRENT_LIB\.samples_file.txt;
fi

### At this point they're ready for any further filtering we want to do.
### For gene expression purposes, we will want to filter out rRNA and ERCC reads,
### But for assembly, we will leave them in.

### Next is library assembly.

mv *.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.err polya_logs/$CURRENT_LIB\_$LIBEXT.polya.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.err
mv *.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.out polya_logs/$CURRENT_LIB\_$LIBEXT.polya.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.out
