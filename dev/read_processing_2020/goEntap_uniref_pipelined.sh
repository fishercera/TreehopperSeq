#!/bin/bash
#SBATCH --job-name=goEntap
#SBATCH --mail-user=cera.fisher@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -c 8 
#SBATCH -N 1
#SBATCH --requeue
#SBATCH -o entap.%N.%A_%a.out
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=24G

source ~/.profile
## Read the libs file to get the array of library prefixes. 
mapfile -t LIB < $1 

CURRENT_LIB=${LIB[$SLURM_ARRAY_TASK_ID]}
LIBEXT=`grep $CURRENT_LIB 2020_library_basenames.txt | cut -f 2`
SPECIES=`grep $CURRENT_LIB 2020_library_information.txt | cut -f 3`
echo $SPECIES
## Need to upload a tab delimited file, 2020_library_info.txt
## that has the necessary info in the correct format 

OUTDIR=$LIBEXT\.entap
#/labs/Wegrzyn/EnTAP/EnTAP_v0.9.0/EnTAP/
module load anaconda/2.4.0
module load diamond/0.9.25
module load interproscan/5.25-64.0
module load rsem/1.2.12
module unload perl
## Check the format for the Alignment bowtie and the Trinity name
#HvitL5Ten_Eye.EJCF002_K11.v1.Trinity.fasta

LIB=$LIBEXT\.$CURRENT_LIB
echo "Processing $LIB"

if [ ! Assemblies/$LIB.v2.Trinity.fasta ]
then 
   echo "Assembly not found." 
else 
/labs/Wegrzyn/EnTAP/EnTAP --runP -d diamondDBs/uniref90.dmnd -i Assemblies/$LIB.v2.Trinity.fasta  -a Alignments/$CURRENT_LIB\.$LIBEXT\.rsem.bowtie2.bam --taxon $SPECIES  -c bacteria -c fungi -c homo_sapiens -c viridiplantae -c archaea --threads 8 --out-dir $OUTDIR  --trim --ontology 0 --ontology 1 --protein pfam
wait
fi 

function EntapFin {
CURRENT_LIB=$1
LIBEXT=$2

  cd $LIBEXT.entap

  sed -i "s/^>TRINITY/>$LIBEXT\_TRINITY/g" final_*notated.f*
  sed -i "s/^TRINITY/$LIBEXT\_TRINITY/g" final_annotations*.tsv
  for file in final_*notated*; do mv $file ${file/final/$LIBEXT}; done
  for file in final_annotations*.tsv; do mv $file ${file/final/$LIBEXT}; done

  cd entap_out

  sed -i "s/^>TRINITY/>$LIBEXT\_TRINITY/g" *.fasta

  cd ..

  mv entap_out $LIBEXT\_entap_out/
  mv $LIBEXT\_*notat* $LIBEXT\_entap_out/.

  cd ..

  touch $CURRENT_LIB\.$LIBEXT\.entap.finalized.ok

}
if [ -s $LIBEXT.entap/final_annotations_lvl4_no_contam.tsv ];
    then
    EntapFin $CURRENT_LIB $LIBEXT
else
    echo "Not finding the annotations file where I expect it" && exit 1
fi

