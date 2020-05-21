#!/bin/bash
#
#SBATCH -p himem,himem2 # partition (queue)
#SBATCH -N 1                      # number of nodes
#SBATCH -J trinity 
#SBATCH -c 16                      # number of cores
#SBATCH --mem=64G
#SBATCH --qos=himem
#SBATCH -o trinity_%N.%A_%a.out        # STDOUT
#SBATCH --mail-type=END,FAIL      # notifications for job done & fail
#SBATCH --mail-user=cera.fisher@uconn.edu # send-to address
source ~/.profile
LIBSLIST=$1
module load bowtie/1.1.2
module load bowtie2/2.2.9
module load samtools/1.9
module load bamtools/2.5.1
module load trinity/2.8.5

## Read the libs file to get the array of library prefixes. 
mapfile -t LIB < $LIBSLIST 

CURRENT_LIB=${LIB[$SLURM_ARRAY_TASK_ID]}
LIBEXT=`grep $CURRENT_LIB 2020_library_basenames.txt | cut -f 2`
#EJCF002_I04.HvitL4.421_Leg.P1.trimmed2.fastq.gz

# The files in SampleFiles properly reference the AssembledReads directory off 
# the /labs/ folder. So no copying of files to my working directory is needed. 

samplesfile=SampleFiles/$CURRENT_LIB.samples_file.txt 
out=$LIBEXT.$CURRENT_LIB.v2


#### Here's where we run Trinity
Trinity --max_memory 64G  --seqType fq --samples_file $samplesfile --CPU 16 --min_contig_length 250 --output $out\.trinity/
####
touch $out\.trinity.finished

cd $out\.trinity/ 

if [ -s Trinity.fasta ]; ## success! 
	then perl /isg/shared/apps/trinity/2.2.0/util/TrinityStats.pl Trinity.fasta > $out\.Trinity.fasta.stats
	mv Trinity.fasta $out\.Trinity.fasta
	mv Trinity.timing $out\.Trinity.timing
	mv Trinity.fasta.gene_trans_map $out\.Trinity.fasta.gene_trans_map
	mv $out\.Trinity.* ../Assemblies/.;
	if [ -s ../Assemblies/$out\.Trinity.fasta ]; then echo "Assembly renamed and moved to Assemblies directory."
		touch ../$out\.safe_to_delete
	fi 
fi 

cd ..

if [ -e $out\.safe_to_delete ]; 
	then rm -rf $out\.trinity/;
fi 

####
# Now we're done, so we can go ahead and move our log files to their directories
mv *.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.err trinity_logs/$LIBEXT.$CURRENT_LIB.trinity.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.err
mv *.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.out trinity_logs/$LIBEXT.$CURRENT_LIB.trinity.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.out


