[![DOI](https://zenodo.org/badge/141181484.svg)](https://zenodo.org/badge/latestdoi/141181484)

# TreehopperSeq
A collection of utility scripts, data analysis pipeline, and documentation for RNA-Seq with non-model organisms.

Software used as part of this repo:
- Trinity (trinityRNAseq) / bowtie / RSEM
- TRIMMOMATIC
- Bowtie2
- EnTAP 
- R (Bioconductor packages)


## RNA-Seq read processing pipeline
**TODO** Add scripts for the RNA QC pipeline and a README. 
1. Quality trim reads _DONE_
2. Strip PolyA tails _DONE_
3. Assemble

Make small test data set

## Annotation and Assembly Refinement
**TODO** Add scripts for using EnTAP to annotate reads, USEARCH to cluster proteomes, and scripts to select nucleotide sequences representative of the clustered proteome (the refined assembly). 

## [GOSeq Walkthrough](https://github.com/fishercera/TreehopperSeq/blob/master/GoSeq_Walkthrough.md)
Files include an R Notebook file outlining the steps used to create the background for GoSeq from EnTAP annotations, create the named vectors that GoSeq will use, and run the enrichment analysis (GoSeq_Walkthrough.RMD). Also included are data files needed to run the script, and the output of running the script with those data files. The file GoTermsMap.py is used to create the one-to-one gene id to GO term map from EnTAP annotations. 

## [TPM Normalization and Between-Species Scaling](https://github.com/fishercera/TreehopperSeq/blob/master/TPM_normalization_scaling.md)
Includes an R Notebook file (should be opened in RStudio) explaining the principles and formula for scaling transcripts per million from one species' transcriptome to another's in order to make the gene expression of the two species more fairly comparable. Also includes **TODO** data files to use to test the code and examine the results of scaling.

## [MLSeq/PLDA classifier Walkthrough](https://github.com/fishercera/TreehopperSeq/blob/master/MLSeq_SpeciesSignal_tuneLength100.md)
This procedure uses scaled TPM of two species for shared single-copy orthologues to find the minimum set of genes that are sufficient to distinguish one species from the other. Files include an R Notebook file with the steps involved and code to examine the results of removing biased genes, a standalone R script meant to be run on a highly parallel compute cluster, and a small set of test data. 
