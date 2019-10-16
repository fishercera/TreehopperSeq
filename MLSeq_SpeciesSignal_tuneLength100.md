Linear discriminant analysis with MLSeq
================
Cera Fisher

``` r
library("ggplot2")
library("DESeq2")
library("MLSeq")
library("dplyr")

CountsData <- read.table("C:/Users/cruth/Google Drive/Treehoppers/ResearchFiles/RNASeq/OrthoFinder/MLSeq/SCOsTPM_2.matrix", header=TRUE, sep="\t")
CountsData_bak <- CountsData
head(CountsData)
colnames(CountsData)
outputPrefix <- "HV_EC_OrthosNEW"
OrthoID.vector <- CountsData$OrthoID
HVid.vector <- CountsData$HVid
ECid.vector <- CountsData$ECid
# It is absolutely critical that the columns of the count matrix and the rows of
# the column data (information about samples) are in the same order. DESeq2 will
# not make guesses as to which column of the count matrix belongs to which row
# of the column data, these must be provided to DESeq2 already in consistent
# order.


##Which Samples Are Wanted##
colData <- read.table("C:/Users/cruth/Google Drive/Treehoppers/ResearchFiles/RNASeq/OrthoFinder/SampleInformation_colData.txt", 
                      sep="\t", header=TRUE) # Read in which samples are wanted
colData$Sample
CountsData <- CountsData[,2:46] # Selecting only the sample information we want. 



## OK do we have all the samples in the right order??? ##
colData$Sample == colnames(CountsData)
## YEP!  They're identical (all "TRUE")

cts <- as.matrix(CountsData)
storage.mode(cts) = "integer"


## Now we can make the DDS object
dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = colData, 
                              design = ~ Tissue)
```

Read in selected markers, filter the counts
matrices.

``` r
selectedMarkers <- read.table("C:/Users/cruth/Google Drive/Treehoppers/ResearchFiles/RNASeq/OrthoFinder/MLSeq/New_SCOs_Tuned_MergedCounts_SelectedMarkers_tuneLength100.txt")

unselected <- read.table("C:/Users/cruth/Google Drive/Treehoppers/ResearchFiles/RNASeq/OrthoFinder/MLSeq/New_SCOs_TunedMergedCounts_UnselectedGenes_tuneLength100.txt")
```

#### Plot the *E. carinata* TPM vs. *H. vitripennis* TPM for genes identified as species markers.

``` r
plotting.data <- data.frame(selectedMarkers[,2:46])
p1 <- ggplot(plotting.data, aes(x=ECty4_Abd, y=HV102_Abd, color="abd")) + 
  geom_point() + 
  geom_point(plotting.data, mapping = aes(x=ECty4_Eye,
                                          y=HV102_Eye,
                                          color="eyes")) +
  geom_point(plotting.data, mapping = aes(x=ECty4_Leg,
                                          y=HV102_Leg,
                                          color="leg")) + 
  geom_point(plotting.data, mapping = aes(x=ECty4_Meso,
                                          y=HV102_Meso,
                                          color="meso")) +
  geom_point(plotting.data, mapping = aes(x=ECty4_Pro,
                                          y=HV102_Pro,
                                          color="pro")) +
  geom_point(plotting.data, mapping = aes(x=ECty4_Ovi,
                                          y=HV102_Ovi,
                                          color="ovi")) + 
  geom_point(plotting.data, mapping = aes(x=ECty4_Wing2,
                                          y=HV102_Wing2,
                                          color="wing2")) +
  geom_point(plotting.data, mapping = aes(x=ECty4_Wing3,
                                          y=HV102_Wing3,
                                          color="wing3")) +
  labs(x="Ecar expression", y="Hvit expression", 
       title="356 transcripts that discriminate for species")


plot <- p1 + coord_cartesian(xlim=c(0, 15000), 
                             ylim=c(0, 15000))
plot
```

![](MLSeq_SpeciesSignal_tuneLength100_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

This plot shows a trident pattern. Extreme outliers in expression along
the axes indicate genes highly expressed in one species and not
expressed in the
other.

#### Plot *E. carinata* TPM vs. *H. vitripennis* TPM of genes not selected as species markers.

``` r
plotting.data2 <- data.frame(unselected)
p2 <- ggplot(plotting.data2, aes(x=ECty4_Abd, y=HV102_Abd,
                                 color="abd")) + 
  geom_point() + 
  geom_point(plotting.data2, mapping = aes(x=ECty4_Eye,
                                           y=HV102_Eye,
                                           color="eyes")) +
  geom_point(plotting.data2, mapping = aes(x=ECty4_Leg,
                                           y=HV102_Leg,
                                           color="leg")) + 
  geom_point(plotting.data2, mapping = aes(x=ECty4_Meso,
                                           y=HV102_Meso,
                                           color="meso")) +
  geom_point(plotting.data2, mapping = aes(x=ECty4_Pro,
                                           y=HV102_Pro,
                                           color="pro")) +
  geom_point(plotting.data2, mapping = aes(x=ECty4_Ovi,
                                           y=HV102_Ovi,
                                           color="ovi")) + 
  geom_point(plotting.data2, mapping = aes(x=ECty4_Wing2,
                                           y=HV102_Wing2,
                                           color="wing2")) +
  geom_point(plotting.data2, mapping = aes(x=ECty4_Wing3,
                                           y=HV102_Wing3,
                                           color="wing3")) +
  labs(x="Ecar expression", y="Hvit expression", 
       title="7279 transcripts classified as unbiased")

plot <- p2 + coord_cartesian(xlim=c(0, 15000), 
                             ylim=c(0, 15000))

plot
```

![](MLSeq_SpeciesSignal_tuneLength100_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

This plot does not have a trident pattern. Most genes with expression
biased towards one species have been removed from the set.

``` r
write.table(selectedMarkers, "MLSeq_PLDA_SelectedMarkers_NEW_356.txt")
write.table(unselected, "UnselectedTranscripts_NEW_7279.txt")
```

### Code used to create PLDA classifier and separate species markers from other genes.

    #~~~~ Cera Fisher (2018) MIT License, use what you like
    ### MLSeq - finding species marker genes with machine learning
    options(stringsAsFactors = FALSE)
    library("dplyr")
    library("DESeq2")
    library("MLSeq")
    
    #Read in TPM matrix
    MergedCounts <- read.delim("SCOsTPM_2.matrix", sep="\t", header=TRUE)
    
    cts <- as.matrix(MergedCounts[,2:46])
    storage.mode(cts) = "integer"
    
    class <- data.frame(condition = factor(rep(c("Ecar","Hvit"),c(24,21))))
    ## Setting up a Class object for DESeq2
    
    
    set.seed(2128)
    vars <- sort(apply(cts, 1, var, na.rm = TRUE), decreasing=TRUE)
    data <- cts # Operating on the whole data set 
    
    ## You can randomly select the set of test samples, but with a small sample set, it's 
    ## possible to get too many of one species; I chose to hand curate the test samples, 
    ## but the below two lines will do it randomly. 
    # nTest <- ceiling(ncol(data) * 0.3)
    # ind <- sample(ncol(data), nTest, FALSE)
    
     GoodInd <- read.table("good_ind2.txt", sep="\t") # using the same set of samples as the old data 
     ind <- as.vector(GoodInd)
     ind <- ind[,1]
    data.train <- as.matrix(data[ ,-ind] + 1) 
    data.test <- as.matrix(data[ ,ind] + 1) 
    classtr <- data.frame(condition = class[-ind, ]) 
    classts <- data.frame(condition = class[ind, ])
    cts.train <- as.matrix(cts[ ,-ind] + 1)
    
    library("DESeq2")
    cts.trainS4 <- DESeqDataSetFromMatrix(countData = cts.train, colData=classtr, design=formula(~condition))
    featureData <- data.frame(gene=MergedCounts$OrthoID)
    
    mcols(cts.trainS4) <- DataFrame(mcols(cts.trainS4), featureData)
    mcols(cts.trainS4) <- DataFrame(mcols(cts.trainS4), data.frame(HVid=MergedCounts$HVid))
    mcols(cts.trainS4) <- DataFrame(mcols(cts.trainS4), data.frame(ECid=MergedCounts$ECid))
    ctrl.PLDA <- discreteControl(method="repeatedcv", number=30, tuneLength=100, repeats=100000, parallel=TRUE)
    
    ## We're setting tuneLength=100 to let the classifier try multiple different tuning parameters
    ## (rho) -- it will settle on the one that gives the sparsest model with the highest accuracy
    
    
    fit.all.PLDA <- classify(cts.trainS4, method="PLDA", preProcess="deseq-vst", 
      control=discreteControl(ctrl.PLDA))
    
    png("SCOs_Filt.TuneLength100.fit.all.PLDA.plot.png")
    plot(fit.all.PLDA)
    dev.off()
    
    
    
    trained(fit.all.PLDA)
    
    ### During my run -- 
    # The optimum model is obtained when rho = 24.08461 with an overall accuracy of
    # Accuracy = 0.9730 over folds. On the average 320.48 out of 7635 features was used
    # in the classifier.
    
    
    Markers <- selectedGenes(fit.all.PLDA)
    Markers.Counts <- MergedCounts[Markers,]
     write.table(Markers.Counts, "New_SCOs_Tuned_MergedCounts_SelectedMarkers_tuneLength100.txt")
     
     UnselectedGenes <- MergedCounts[(which(!(MergedCounts$OrthoID %in% MergedCounts[Markers, ]$OrthoID))), ]
     write.table(UnselectedGenes, "New_SCOs_TunedMergedCounts_UnselectedGenes_tuneLength100.txt", sep="\t", quote=FALSE)
