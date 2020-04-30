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

2)  Calculate &#945;

<!-- end list -->

``` r
HvitN1 <- as.numeric(length(Hvit.TPM$HVid))
                                            # 19,126
EcarN2 <- as.numeric(length(Ecar.TPM$ECid))
                                            # 19,975

sumN2.EcarTPM <- (sum(Ecar.TPM[,2:42]))/41
sumN2.EcarTPM
```

    ## [1] 1e+06

``` r

## Sum of TPM for any given sample should, by definition, be 1,000,000
# The average TPM for Ecar is the sum divided by the number of transcripts. 
Ec.avg.TPM <- sumN2.EcarTPM/EcarN2



## Getting the sum of Ecar TPM for the # of transcripts in Hvit's transcriptome
## i.e, multiplying the average times 19,126
sum.ECavgTPM.HvitN1 <- Ec.avg.TPM * HvitN1

# To get alpha, divide that amount by 1,000,000
alpha <- sum.ECavgTPM.HvitN1 * (10**(-6))

alpha
```

    ## [1] 0.9574969

``` r
beta <- sum.ECavgTPM.OfasN3 * (10**(-6))
```

This value for 	&#945;, 0.957….etc is very close to the ratio of the
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

Multiplying by 	&#945; results in our Hvit TPM numbers being just a
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

Now, let’s save our scaled and size-factor normalized TPM to new files
to use later.

``` r
write.table(Hvit.TPM.scaled, "Hvit.90.isoforms.TPM.scaled.matrix", sep="\t")

library(DESeq2)
```
We loaded DESeq2 in order to normalize the scaled values by library size within *H. vitripennis*
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
```

And then we'll want to do the same for *E. carinata*

``` r
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
