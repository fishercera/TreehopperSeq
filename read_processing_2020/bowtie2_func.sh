#!/usr/bin/bash

##############################################################################
# This script runs libraries through the bowtie2 decontamination routine.    #
# USAGE: bowtie2_func.sh <lib-base-name> <List that is the names of indices> #
# For example: bowtie2_func.sh ELJ1393 Greg_rRNA dmelan scerev               #
##############################################################################

source ~/.profile                                  # My bash .profile provides many useful environment variables.

# filterContamsPE - this works on paired reads and filters out the pairs that align concordantly. We're being very generous
# with the definition of concordant, allowing dovetailed reads just in case we get readthrough.
#USAGE: filterContamsPE <lib-base-name> <index> <index-path>

function filterContamsPE {
lib=$1 # The library base, ie ELJ1393
index=$2 # The index base, ie Greg_rRNA
#indexPath=$3 # I have bowtie2 indices in the local directory bt2, but you could put them somewhere else
              # and uncomment this variable, and then pass the directory info along to the function.
indexPath=bt2

echo "Starting bowtie2 with parameters $lib, $index, $indexPath"
unconc=input/scratch/$lib\.$index.P%.filtered.fastq.gz
alconc=input/scratch/$lib\.$index.P%.contams.fastq.gz

bowtie2 -q --phred33  -p 4 --very-fast  --dovetail --al-conc-gz $alconc --un-conc-gz $unconc -x $indexPath/$index -1 <(zcat input/$lib.R1.fastq.gz) -2 <(zcat input/$lib.R2.fastq.gz) -S output/$lib\.$index\.bam


# These lines will ensure that the next run through this is going to use the filtered files, and we'll keep filtering stuff out.
      cp input/scratch/$lib\.$index.P1.filtered.fastq.gz input/$lib.R1.fastq.gz
      cp input/scratch/$lib\.$index.P2.filtered.fastq.gz input/$lib.R2.fastq.gz
      touch $lib\.Pairedbowtie$index.done

}


function filterContamsSE {
lib=$1 # The library base, ie ELJ1393
index=$2 # The index base, ie Greg_rRNA
##indexPath=$3 # The path to the index, ie ~/mini/bt2
indexPath=bt2

echo "Starting bowtie2 with parameters $lib, $index, $indexPath"
un=input/scratch/$lib\.$index\.U1.filtered.fastq.gz
al=input/scratch/$lib\.$index\.U1.contams.fastq.gz
##unconc=input/scratch/$lib\.$index.P%.filtered.fastq.gz
##alconc=input/scratch/$lib\.$index.P%.contams.fastq.gz
bowtie2 -q --phred33  -p 4 --very-fast  --al-gz $al --un-gz $un -x $indexPath/$index -U <(zcat input/$lib.U1.fastq.gz) \
 -S output/$lib\.$index\.U1.bam
wait;

un=input/scratch/$lib\.$index\.U2.filtered.fastq.gz
al=input/scratch/$lib\.$index\.U2.contams.fastq.gz
bowtie2 -q --phred33 -p 4 --very-fast --al-gz $al --un-gz $un -x $indexPath/$index -U <(zcat input/$lib.U2.fastq.gz) \
 -S output/$lib\.$index\.U2.bam

# These lines will ensure that the next run through this is going to use the filtered files, and we'll keep filtering stuff out.
      cp input/scratch/$lib\.$index\.U1.filtered.fastq.gz input/$lib.U1.fastq.gz
      cp input/scratch/$lib\.$index\.U2.filtered.fastq.gz input/$lib.U2.fastq.gz
      touch $lib\.SEbowtie.$index.done

}


args="$@"
numargs=$#
echo $args
echo $numargs
echo "libbase = $1"
lib=$1
indexPath=bt2
shift
echo "indices = $@"
indices=$@
echo $indices

module load bowtie2

for i in $indices

do
  echo $i
  echo "--filtering--"
  filterContamsPE $lib $i $indexPath
  PECOUNT=`zgrep -c '^+$' input/$lib.R1.fastq.gz`
  filterContamsSE $lib $i $indexPath
  U1COUNT=`zgrep -c '^+$' input/$lib.U1.fastq.gz`
  U2COUNT=`zgrep -c '^+$' input/$lib.U2.fastq.gz`
  echo "Remaining PE reads: $PECOUNT"
  echo "Remaining U1 reads: $U1COUNT"
  echo "Remaining U2 reads: $U2COUNT"
  echo "--filtered--"
  echo $lib
done


#cp input/scratch/$lib.U.filtered.fastq output/.
cp input/$lib.R1.fast* output/.
cp input/$lib.R2.fast* output/.
cp input/$lib.U*.fastq.gz output/.



## Removing the BAM files because it's too large.
#EJCF002_N01.EbinL5T_Leg.EbinrRNA.bam
rm -f output/$lib\*.bam
##### Move the contamreads files to ContamReads
mv input/scratch/$lib\.*.contams.fastq.gz $LABS_DIR/ContamReads/.
##### Move filtered reads to FilteredReads/; these are for gene expression analysis
if [ -s output/$lib\.R1.fastq.gz ]; then
    mv output/$lib\.*.fastq.gz $LABS_DIR/FilteredReads/.
else
   echo "I'm not finding the expected files in output/ ." && exit 1
fi
echo "Contaminant alignments in BAM files at $LABS_DIR/ContamReads/."
echo "Filtered reads at $LABS_DIR/FilteredReads/."