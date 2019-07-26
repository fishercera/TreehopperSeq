TPM Normalization Between Species
================
Cera Fisher
September 14, 2018

## Scaling TPM between two species with different numbers of annotated transcripts

*based on Musser & Wagner (2015), JEZ:B*

TPM, transcripts per million mapped, is a way of normalizing RNAseq data
that accounts for differences in library size by scaling the abundance
of a transcript (the “counts”) to the total number of transcripts
assumed to be present in the transcriptome. In short, TPM = count for
transcript “i” / total number of annotated transcripts \* a scaling
factor.

Necessarily, TPM from one species does not map to the TPM of another
species if their transcriptomes are of different sizes, which they
almost certainly are, so the values for the species with the smaller
transcriptome will be inflated relative to the species with the larger
transcriptome. According to Musser & Wagner, these values can be
rescaled by calculating a scaling factor &#945;:



![This would all look quite nice when rendered to HTML. For github, 
we're stuck with an image of the formula.](https://i.imgur.com/yDpZjwm.png)

Where N1 = the number of transcripts for species B (the smaller set),
*j* = all the transcripts in the set, and tpm(Aj) is the transcripts per
million for species A for each of the genes *j* in the set.

Here is how I interpret this to work in my transcriptomes for *Entylia
carinata* and *Homalodisca vitripennis*.

1)  Read in the un-normalized TPM values from RSEM/edgeR.

<!-- end list -->

``` r
library("dplyr")
options(stringsAsFactors = FALSE)
Ecar.TPM <- read.table("C:/Users/cruth/Treehoppers/ResearchFiles/RNASeq/GeneExpression_2018/Annotation/ECEF_Refined.isoforms.TPM.not_cross_norm", sep="\t", header=TRUE)
colnames(Ecar.TPM)
```

    ##  [1] "X"            "ECA_Abd"      "ECA_Ovi"      "ECB_Abd"     
    ##  [5] "ECC_Abd"      "ECC_Ovi"      "ECEF_Abd"     "ECEF_Eye"    
    ##  [9] "ECEF_Leg"     "ECEF_Meso"    "ECEF_Ovi"     "ECEF_Pro"    
    ## [13] "ECEF_Wing2"   "ECEF_Wing3"   "ECFisC_Abd"   "ECFisC_Eye"  
    ## [17] "ECFisC_Leg"   "ECFisC_Meso"  "ECFisC_Ovi"   "ECFisC_Pro"  
    ## [21] "ECFisC_Wing2" "ECFisC_Wing3" "ECLegs_A"     "ECLegs_B"    
    ## [25] "ECLegs_C"     "ECMeso_A"     "ECMeso_B"     "ECMeso_C"    
    ## [29] "ECPro_A"      "ECPro_B"      "ECPro_C"      "ECty4_Abd"   
    ## [33] "ECty4_Eye"    "ECty4_Leg"    "ECty4_Meso"   "ECty4_Ovi"   
    ## [37] "ECty4_Pro"    "ECty4_Wing2"  "ECty4_Wing3"  "ECWings_A"   
    ## [41] "ECWings_B"    "ECWings_C"

``` r
colnames(Ecar.TPM)[1] <- "ECid"

Hvit.TPM <- read.table("C:/Users/cruth/Treehoppers/ResearchFiles/RNASeq/GeneExpression_2018/Annotation/HV_Refined.isoform.TPM.not_cross_norm", sep="\t", header=TRUE)
colnames(Hvit.TPM)
```

    ##  [1] "X"           "HV102_Abd"   "HV102_Eye"   "HV102_Leg"   "HV102_Meso" 
    ##  [6] "HV102_Ovi"   "HV102_Pro"   "HV102_Wing2" "HV102_Wing3" "HV3a_Abd"   
    ## [11] "HV3a_Eye"    "HV3a_Leg"    "HV3a_Meso"   "HV3a_Ovi"    "HV3a_Pro"   
    ## [16] "HV3a_Wing2"  "HV3a_Wing3"  "HV7_Abd"     "HV7_Leg"     "HV7_Meso"   
    ## [21] "HV7_Pro"     "HV7_Wing2"

``` r
colnames(Hvit.TPM)[1] <- "HVid"

Ofas.TPM <- read.table("C:/Users/cruth/Treehoppers/ResearchFiles/SRAProject/Ofasciatus/Ofas90.isoform.TMM.EXPR.matrix", sep="\t", header=TRUE)

colnames(Ofas.TPM)
```

    ##  [1] "X"                 "Ofas_T1S1"         "Ofas_T1S2"        
    ##  [4] "Ofas_T1S3"         "Ofas_T1SCR_RNAi_1" "Ofas_T1SCR_RNAi_2"
    ##  [7] "Ofas_T2S1"         "Ofas_T2S2"         "Ofas_T2S3"        
    ## [10] "Ofas_T3S1"         "Ofas_T3S2"         "Ofas_T3S3"

``` r
colnames(Ofas.TPM)[1] <- "OFid"
```

2)  Calculate \(\alpha\)

<!-- end list -->

``` r
HvitN1 <- as.numeric(length(Hvit.TPM$HVid))
                                            # 19,126
EcarN2 <- as.numeric(length(Ecar.TPM$ECid))
                                            # 19,975
OfasN3 <- as.numeric(length(Ofas.TPM$OFid))
                                            # 19,811

sumN2.EcarTPM <- (sum(Ecar.TPM[,2:42]))/41
sumN2.EcarTPM
```

    ## [1] 1e+06

``` r
sumN3.OfasTPM <- sum(Ofas.TPM[,3]) ## Ofas T1S1 is inflated. :( 

sumN3.OfasTPM <- sum(Ofas.TPM[,2:12])/11 ## Ofas T1S1 is inflated. :( 

## Sum of TPM for any given sample should, by definition, be 1,000,000
# The average TPM for Ecar is the sum divided by the number of transcripts. 
Ec.avg.TPM <- sumN2.EcarTPM/EcarN2



## Getting the sum of Ecar TPM for the # of transcripts in Hvit's transcriptome
## i.e, multiplying the average times 19,126
sum.ECavgTPM.HvitN1 <- Ec.avg.TPM * HvitN1
sum.ECavgTPM.OfasN3 <- Ec.avg.TPM * OfasN3

# To get alpha, divide that amount by 1,000,000
alpha <- sum.ECavgTPM.HvitN1 * (10**(-6))

alpha
```

    ## [1] 0.9574969

``` r
beta <- sum.ECavgTPM.OfasN3 * (10**(-6))
```

This value for \(\alpha\), 0.957….etc is very close to the ratio of the
smaller number of transcripts to the larger:

``` r
checksum <- HvitN1/EcarN2
checksum
```

    ## [1] 0.9574969

``` r
checksum - alpha
```

    ## [1] 7.006075e-09

Which perhaps should be expected, since those are the only values in
this calculation that don’t cancel out.

Scaling Hvit.TPM, then, goes like this

``` r
HV_scaled <- Hvit.TPM[,-1]
row.names(HV_scaled) <- Hvit.TPM[,1]
Hvit.TPM.scaled <- HV_scaled * alpha
dim(HV_scaled)                       
```

    ## [1] 19126    21

``` r
head(Hvit.TPM[,3])
```

    ## [1] 12.68  3.08 30.92  2.55  0.29  6.34

``` r
head(Hvit.TPM.scaled[,2])
```

    ## [1] 12.1410602  2.9490903 29.6058030  2.4416170  0.2776741  6.0705301

``` r
OF_scaled <- Ofas.TPM[,-1]
row.names(OF_scaled) <- Ofas.TPM[,1]
Ofas.TPM.scaled <- OF_scaled * beta
dim(OF_scaled)                       
```

    ## [1] 15627    11

``` r
head(Ofas.TPM[,3])
```

    ## [1]  64.857   0.000   1.572   6.955   0.539 104.392

``` r
head(Ofas.TPM.scaled[,2])
```

    ## [1] 50.7394409  0.0000000  1.2298195  5.4410906  0.4216747 81.6687746

Multiplying by \(\alpha\) results in our Hvit TPM numbers being just a
little smaller, though it will matter a lot for some of the outrageously
large numbers–

``` r
which(Hvit.TPM[,3] == max(Hvit.TPM[,3]))
```

    ## [1] 356

``` r
#356
Hvit.TPM[,3][356]
```

    ## [1] 136409.6

``` r
Hvit.TPM.scaled[,2][356]
```

    ## [1] 130611.7

``` r
which(Ofas.TPM[,3] == max(Ofas.TPM[,3]))
```

    ## [1] 5993

``` r
Ofas.TPM[,3][5993]
```

    ## [1] 39721.33

``` r
Ofas.TPM.scaled[,2][5993]
```

    ## [1] 31075.11

Now, let’s save our scaled and size-factor normalized TPM to new files
to use later.

``` r
head(Ofas.TPM.scaled)
```

    ##                                             Ofas_T1S1  Ofas_T1S2
    ## Ofas_T3S2_TRINITY_DN14866_c4_g5_i1          0.2151402 50.7394409
    ## Ofas_T3S2_TRINITY_DN19291_c0_g1_i1          0.0000000  0.0000000
    ## Ofas_T3S3_TRINITY_DN10610_c0_g2_i1          0.0000000  1.2298195
    ## Ofas_T1SCR_RNAi_1_TRINITY_DN22761_c0_g1_i1  1.5466623  5.4410906
    ## Ofas_T3S2_TRINITY_DN15939_c0_g1_i1          1.8048305  0.4216747
    ## Ofas_T1S1_TRINITY_DN5502_c0_g1_i1          52.8884956 81.6687746
    ##                                             Ofas_T1S3 Ofas_T1SCR_RNAi_1
    ## Ofas_T3S2_TRINITY_DN14866_c4_g5_i1         47.5350258         6.9275136
    ## Ofas_T3S2_TRINITY_DN19291_c0_g1_i1          0.5851813         0.1901057
    ## Ofas_T3S3_TRINITY_DN10610_c0_g2_i1          1.6397593         1.0678776
    ## Ofas_T1SCR_RNAi_1_TRINITY_DN22761_c0_g1_i1  7.6949773         6.9768002
    ## Ofas_T3S2_TRINITY_DN15939_c0_g1_i1          1.8822809         0.7369529
    ## Ofas_T1S1_TRINITY_DN5502_c0_g1_i1          96.7098109       119.8189948
    ##                                            Ofas_T1SCR_RNAi_2   Ofas_T2S1
    ## Ofas_T3S2_TRINITY_DN14866_c4_g5_i1                  2.943900 816.8622048
    ## Ofas_T3S2_TRINITY_DN19291_c0_g1_i1                  2.048917   0.9348818
    ## Ofas_T3S3_TRINITY_DN10610_c0_g2_i1                  0.000000   0.0000000
    ## Ofas_T1SCR_RNAi_1_TRINITY_DN22761_c0_g1_i1          1.667923   0.7080068
    ## Ofas_T3S2_TRINITY_DN15939_c0_g1_i1                  2.548042   0.8081447
    ## Ofas_T1S1_TRINITY_DN5502_c0_g1_i1                  59.788628  53.3500691
    ##                                             Ofas_T2S2   Ofas_T2S3
    ## Ofas_T3S2_TRINITY_DN14866_c4_g5_i1           1.107776   0.4115045
    ## Ofas_T3S2_TRINITY_DN19291_c0_g1_i1           1.456695   1.6953046
    ## Ofas_T3S3_TRINITY_DN10610_c0_g2_i1           0.000000   0.0000000
    ## Ofas_T1SCR_RNAi_1_TRINITY_DN22761_c0_g1_i1   2.442428   0.0000000
    ## Ofas_T3S2_TRINITY_DN15939_c0_g1_i1           1.884628   3.1879862
    ## Ofas_T1S1_TRINITY_DN5502_c0_g1_i1          149.752425 146.8859756
    ##                                              Ofas_T3S1  Ofas_T3S2
    ## Ofas_T3S2_TRINITY_DN14866_c4_g5_i1         295.4516093 69.2360195
    ## Ofas_T3S2_TRINITY_DN19291_c0_g1_i1           1.1476750  7.1629943
    ## Ofas_T3S3_TRINITY_DN10610_c0_g2_i1           1.9182680  0.1439483
    ## Ofas_T1SCR_RNAi_1_TRINITY_DN22761_c0_g1_i1   0.8652547  0.0000000
    ## Ofas_T3S2_TRINITY_DN15939_c0_g1_i1           1.3236988  2.9744107
    ## Ofas_T1S1_TRINITY_DN5502_c0_g1_i1           77.9980920 62.4915707
    ##                                             Ofas_T3S3
    ## Ofas_T3S2_TRINITY_DN14866_c4_g5_i1         156.244182
    ## Ofas_T3S2_TRINITY_DN19291_c0_g1_i1           5.597556
    ## Ofas_T3S3_TRINITY_DN10610_c0_g2_i1           9.110991
    ## Ofas_T1SCR_RNAi_1_TRINITY_DN22761_c0_g1_i1   0.000000
    ## Ofas_T3S2_TRINITY_DN15939_c0_g1_i1           1.450436
    ## Ofas_T1S1_TRINITY_DN5502_c0_g1_i1           27.710837

``` r
write.table(Ofas.TPM.scaled, "Ofas.90.isoforms.TPM.scaled.matrix", sep="\t")

write.table(Hvit.TPM.scaled, "Hvit.90.isoforms.TPM.scaled.matrix", sep="\t")

library(DESeq2)
```

    ## Warning: package 'DESeq2' was built under R version 3.5.2

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colMeans, colnames, colSums, dirname, do.call, duplicated,
    ##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
    ##     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
    ##     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
    ##     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which, which.max,
    ##     which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Warning: package 'GenomeInfoDb' was built under R version 3.5.2

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## Warning: package 'matrixStats' was built under R version 3.5.3

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## Loading required package: BiocParallel

    ## Warning: package 'BiocParallel' was built under R version 3.5.2

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply

``` r
colData <- read.table("SampleInformation_colData.txt", header=TRUE)
hvCol <- colData[25:45,]
hvTPM <- as.matrix(Hvit.TPM.scaled)
storage.mode(hvTPM) = "integer"
hv.dds <- DESeqDataSetFromMatrix(hvTPM, colData = hvCol, 
                                 design = ~ Tissue + Pool)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
hv.dds <- estimateSizeFactors(hv.dds)
hv.dds <- estimateDispersions(hv.dds)
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
hvScaledNorm <- as.data.frame(assay(hv.dds), normalized = TRUE)
write.table(hvScaledNorm, "Hvit_Scaled_SizeNormed_Integer_TPM.txt", sep="\t")

ECid <- Ecar.TPM[,1]
Ecar.TPM <- Ecar.TPM[,-1]
Ecar.TPM <- Ecar.TPM[,c(6:13,14:21,31:38)]
rownames(Ecar.TPM) <- ECid
ec.dds <- as.matrix(Ecar.TPM)
ecCol <- colData[1:24,]
storage.mode(ec.dds) = "integer"

ec.dds <- DESeqDataSetFromMatrix(ec.dds, colData = ecCol,
                                 design = ~ Tissue + Pool)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
ec.dds <- estimateSizeFactors(ec.dds)
ec.dds <- estimateDispersions(ec.dds)
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## -- note: fitType='parametric', but the dispersion trend was not well captured by the
    ##    function: y = a/x + b, and a local regression fit was automatically substituted.
    ##    specify fitType='local' or 'mean' to avoid this message next time.

    ## final dispersion estimates

``` r
ecNorm <- as.data.frame(assay(ec.dds), normalized=TRUE)
write.table(ecNorm, "Ecar_SizeNormed_Integer_TPM.txt", sep="\t")
```
