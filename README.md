# TreehopperSeq
A collection of utility scripts, data analysis pipeline, and documentation for RNA-Seq with non-model organisms.

Software used as part of this repo:
- Trinity (trinityRNAseq)
- TRIMMOMATIC
- EnTAP 
- R (Bioconductor packages)


## GOSeq Walkthrough
Files include an R Notebook file outlining the steps used to create the background for GoSeq from EnTAP annotations, create the named vectors that GoSeq will use, and run the enrichment analysis (GoSeq_Walkthrough.RMD). Also included are data files needed to run the script, and the output of running the script with those data files. The file GoTermsMap.py is used to create the one-to-one gene id to GO term map from EnTAP annotations. 

