#!/usr/bin/bash
#SBATCH -p general # partition (queue)
#SBATCH -N 1                      # number of nodes
#SBATCH -J RSEM
#SBATCH -n 4                      # number of cores
#SBATCH --mem=8G
#SBATCH --qos=general
#SBATCH -o RSEM.%N.%A_%a.out        # STDOUT
#SBATCH -e RSEMstats.%N.%A_%a.err
#SBATCH --mail-type=END,FAIL      # notifications for job done & fail
#SBATCH --mail-user=cera.fisher@uconn.edu # send-to address

############################# Set up environment #################################

source ~/.profile #acquire special environmental variables [$LABS_DIR, $LINUXSHARE]
## Read the libs file to get the array of library prefixes. 
mapfile -t LIB < $1

CURRENT_LIB=${LIB[$SLURM_ARRAY_TASK_ID]}
LIBEXT=`grep $CURRENT_LIB ../2020_library_basenames.txt | cut -f 2`
RNADB=`grep $CURRENT_LIB ../2020_library_information.txt | cut -f 6`
TRANSCRIPTS=$LIBEXT\.$CURRENT_LIB\.v2.Trinity.fasta
WORKINGDIR=~/bio/Membracids/ReadProcessing/
module load bowtie2/2.3.5.1
module load rsem/1.3.0
module load trinity/2.8.5
module load samtools/1.9
echo "$TRINITY_HOME"


########################## Script Follows ######################################

if [ ! -f $LABS_DIR/PreEntap_Assemblies/$TRANSCRIPTS ]; 
then 
    echo "Cannot find $TRANSCRIPTS in $LABS_DIR/PreEntap_Assemblies/, exiting." && exit 1
else 
    cp $LABS_DIR/PreEntap_Assemblies/$TRANSCRIPTS .
fi

if [ ! -s $LABS_DIR/FilteredReads/$CURRENT_LIB\.$LIBEXT.R1.fastq.gz ];
then
    echo "Read files are not where I expect in $LABS_DIR/FilteredReads." && exit 1
fi

LEFT=$LABS_DIR/FilteredReads/$CURRENT_LIB\.$LIBEXT\.R1.fastq.gz 
RIGHT=$LABS_DIR/FilteredReads/$CURRENT_LIB\.$LIBEXT\.R2.fastq.gz

#### Running the Trinity-included align_and_estimate_abundance script
#### We need to end up with an RSEM valid BAM file, 
#### And rather than reinvent the wheel, 
#### Let's just use the script Brian Haas wrote that we know works. 

perl $TRINITY_HOME/util/align_and_estimate_abundance.pl \
 --transcripts $TRANSCRIPTS --est_method RSEM --output_dir $LIBEXT\.out --aln_method bowtie2 \
 --left $LEFT --right $RIGHT \
--seqType fq --thread_count 4 --trinity_mode --prep_reference 

cd $LIBEXT\.out
if [ -f *.bam.ok ]; then
    mv bowtie2.bam $CURRENT_LIB\.$LIBEXT\.bowtie2.bam
    mv bowtie2.bam.for_rsem.bam ~/bio/Membracids/entap/Alignments/$CURRENT_LIB\.$LIBEXT\.rsem.bowtie2.bam
fi

if [ -f RSEM.isoforms.results.ok ]; then
   echo "Results are acquired, let's clean up."
   rm -f ../$TRANSCRIPTS*
   mv RSEM.genes.results $CURRENT_LIB\.$LIBEXT\.genes.results
   mv RSEM.isoforms.results $CURRENT_LIB\.$LIBEXT\.genes.results
fi


mv RSEM.*.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID\.err RSEM_stderr/.
mv RSEMstats.*.$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID\.out RSEM_logs/.



######## trash ############

#if [ ! -s $TRANSCRIPTS.1.bt2 ]; 
#then 
#    PREP_BT2 $TRANSCRIPTS
#fi

#function PREP_BT2 {
#TRANSCRIPTS=$1
#
#bowtie2-build $LABS_DIR/PreEntap_Assemblies/$TRANSCRIPTS $TRANSCRIPTS
#
#}

#function AlignReadsSE {
#lib=$1 # The library basename
#TRANSCRIPTS=$2 # the trinity assembly
#
#if [ ! -d $lib.output ];
#then
#  mkdir $lib.output
#fi
#
#echo "Starting bowtie2 with parameters $lib, $index, $indexPath"
### adding trinity options
#bowtie2 -q --phred33  -p 8 --no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 -x $TRANSCRIPTS -U <(zcat $WORKINGDIR/input/$lib.U.fastq.gz) -S $lib.output/$lib.aln.bam
#
#
#}

