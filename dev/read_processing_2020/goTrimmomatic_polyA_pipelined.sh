#!/bin/bash
#
#SBATCH -p general # partition (queue)
#SBATCH -N 1                      # number of nodes
#SBATCH -J runTrim
#SBATCH -n 8                      # number of cores
#SBATCH --array=[0-215]%8
#SBATCH --mem=16G
#SBATCH --qos=general
#SBATCH -o trimmomatic.%N.%A_%a.out        # STDOUT
#SBATCH --mail-type=END,FAIL      # notifications for job done & fail
#SBATCH --mail-user=cera.fisher@uconn.edu # send-to address
LIBFILE=$1
source ~/.profile
## Make these directories in your workdir before running:
# polyA_output/
# polya_logs/

mapfile -t LIB < $LIBFILE 

CURRENT_LIB=${LIB[$SLURM_ARRAY_TASK_ID]}
LIBEXT=`grep "$CURRENT_LIB" 2020_library_basenames.txt | cut -f 2`

READDIR=output/

wheretrim=/isg/shared/apps/Trimmomatic/0.36
LEFT=$CURRENT_LIB.P1.trimmed.fastq.gz
RIGHT=$CURRENT_LIB.P2.trimmed.fastq.gz
#cp EJCF002_N11_S221_L003_R1_001.fastq.gz ~/bio/Membracids/ReadProcessing/input/.

if [ -s $READDIR/$LEFT ]
then
    module load Trimmomatic
else
   if [ -s $LINUXSHARE/TrimmedReads/$LEFT ]; then 
   cp $LINUXSHARE/TrimmedReads/$LEFT output/$LEFT;
   wait
   cp $LINUXSHARE/TrimmedReads/$RIGHT output/$RIGHT;
   wait
   module load Trimmomatic   
   fi
fi
wait

if [ -s $READDIR/$LEFT ]; then 
    java -jar $wheretrim/trimmomatic-0.36.jar PE -threads 8 $READDIR/$LEFT $READDIR/$RIGHT polyA_output/$CURRENT_LIB.$LIBEXT.P1.trimmed2.fastq.gz polyA_output/$CURRENT_LIB.$LIBEXT.U1.trimmed2.fastq.gz polyA_output/$CURRENT_LIB.$LIBEXT.P2.trimmed2.fastq.gz polyA_output/$CURRENT_LIB.$LIBEXT.U2.trimmed2.fastq.gz ILLUMINACLIP:adapters/polyAtails.fa:2:30:10;
### These will be trimmed2 because I can't think of a better less cumbersome name to mark what part
### of the process they've gone through
else
   echo "Look, the files just aren't here, okay?"
fi

if [ -s polyA_output/$CURRENT_LIB.$LIBEXT.P1.trimmed2.fastq.gz ];
then
   touch $CURRENT_LIB.polyA.done;
   echo "$CURRENT_LIB   $CURRENT_LIB    polyA_output/$CURRENT_LIB.$LIBEXT.P1.trimmed2.fastq.gz  polyA_output/$CURRENT_LIB.$LIBEXT.P2.trimmed2.fastq.gz" > $CURRENT_LIB\.samples_file.txt;
   #mv $READDIR/*$CURRENT_LIB*.fastq.gz /linuxshare/projects/jockuschlab/cfisher/TrimmedReads/.
fi

### At this point they're ready for any further filtering we want to do.
### For gene expression purposes, we will want to filter out rRNA and ERCC reads,
### But for assembly, we will leave them in.

### Next is library assembly.

#mv *.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.err polya_logs/$CURRENT_LIB\_$LIBEXT.polya.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.err
mv *.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.out polya_logs/$CURRENT_LIB\_$LIBEXT.polya.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.out

