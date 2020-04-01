#!/bin/bash
#
#SBATCH -p general # partition (queue)
#SBATCH -N 1                      # number of nodes
#SBATCH -J runTrim
#SBATCH -n 4                      # number of cores
#SBATCH --mem=4G
#SBATCH --qos=general
#SBATCH -o trimmomatic.%N.%A_%a.out        # STDOUT
#SBATCH -e trimmomatic.%N.%A_%a.err        # STDERR
#SBATCH --mail-type=END,FAIL      # notifications for job done & fail
#SBATCH --mail-user=cera.fisher@uconn.edu # send-to address


##### Set up variables
READDIR=/linuxshare/projects/jockuschlab/cfisher/RAW_NovaSeqReads/

mapfile -t LIB < libs.txt

CURRENT_LIB=${LIB[$SLURM_ARRAY_TASK_ID]}
LIBEXT=`grep $CURRENT_LIB 2020_library_basenames.txt | cut -f 2`
 wheretrim=/isg/shared/apps/Trimmomatic/0.36
 LEFT=$CURRENT_LIB\_*_R1_001.fastq.gz
 RIGHT=$CURRENT_LIB\_*_R2_001.fastq.gz
 #cp EJCF002_N11_S221_L003_R1_001.fastq.gz ~/bio/Membracids/ReadProcessing/input/.

 if [ -s $READDIR/$LEFT ];  # if the file exists and has non-zero size in the raw-reads directory
 then                       # copy it to my input directory
     cp $READDIR/$LEFT input/.
     cp $READDIR/$RIGHT input/.
     echo "Raw reads moved from /linuxshare"  # And let me know you did it.
 else                                         # Otherwise,
     echo "$LEFT not found in $READDIR"       # Let me know you couldn't find it.
 fi


# ###### Trimmomatic command commented out for testing purposes, echoing command to file.
 module load Trimmomatic
 java -jar $wheretrim/trimmomatic-0.36.jar PE -threads 8  input/$LEFT input/$RIGHT input/scratch/$CURRENT_LIB.P1.trimmed.fastq.gz input/scratch/$CURRENT_LIB.U1.fastq.gz input/scratch/$CURRENT_LIB.P2.trimmed.fastq.gz input/scratch/$CURRENT_LIB.U2.fastq.gz ILLUMINACLIP:adapters/trimadapters.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:40

 echo "java -jar $wheretrim/trimmomatic-0.36.jar PE -threads 1  input/$LEFT input/$RIGHT input/scratch/$CURRENT_LIB.P1.trimmed.fastq.gz input/scratch/$CURRENT_LIB.U1.fastq.gz input/scratch/$CURRENT_LIB.P2.trimmed.fastq.gz input/scratch/$CURRENT_LIB.U2.fastq.gz ILLUMINACLIP:adapters/trimadapters.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:40"

 if [ -s  input/scratch/$CURRENT_LIB.P1.trimmed.fastq.gz ];
 then
     cp input/scratch/$CURRENT_LIB.P1.trimmed.fastq.gz output/$CURRENT_LIB.P1.trimmed.fastq.gz
     cp input/scratch/$CURRENT_LIB.P2.trimmed.fastq.gz output/$CURRENT_LIB.P2.trimmed.fastq.gz
     cp input/scratch/$CURRENT_LIB.U*.fastq.gz output/.
     touch $CURRENT_LIB\.trimmedreads_moved.ok
     echo "Trimming done and reads moved."
 fi

 echo "$CURRENT_LIB      $CURRENT_LIB    output/$CURRENT_LIB.P1.trimmed.fastq.gz output/$CURRENT_LIB.P2.trimmed.fastq.gz"
 echo "$CURRENT_LIB      $CURRENT_LIB    output/$CURRENT_LIB.P1.trimmed.fastq.gz output/$CURRENT_LIB.P2.trimmed.fastq.gz" > $CURRENT_LIB\.samples_file.txt

 cat $CURRENT_LIB\.samples_file.txt

 mv trimmomatic.*.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.err trimmomatic_logs/$CURRENT_LIB.$LIBEXT.qualtrim.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.err
 mv trimmomatic.*.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.out trimmomatic_logs/$CURRENT_LIB.$LIBEXT.qualtrim.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.out

