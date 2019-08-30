### Cera Fisher 2019
### This code makes use of the heatmap.3 function by Obi Griffith (https://github.com/obigriffith/biostar-tutorials)
### as modified by Brian Haas of the Broad Institute for the Trinity software suite (https://github.com/trinityrnaseq/trinityrnaseq/blob/devel/Analysis/DifferentialExpression/R/heatmap.3.R). 
### The original code was made available under the GPL, and the Trinity software suite is licensed by 
### the Broad Institute under the license available here: https://github.com/trinityrnaseq/trinityrnaseq/blob/devel/LICENSE.txt
### The code below is my own, and licensed under the MIT License. Use what you like. 

library(DESeq2)
library(dplyr)
options(stringsAsFactors = FALSE)
HVdata <- read.table("HV102.genes.counts.matrix.collapsed.wing.isos.txt", sep="\t", header=T, row.names = "X")

HVdata$id <- rownames(HVdata)
col_ordering = c(1,9,17,
                 2,10,
                 3,11,18,
                 4,12,19,
                 5,13,
                 6,14,20,
                 7,15,21,
                 8,16)

Hvit.gene.orthos <- read.table("CandidateWingRelatedGenes_OrthoTranscriptIDs.txt", sep="\t", header=TRUE)
Hvit.final.id <- gsub("_i\\d+", "", Hvit.gene.orthos$Homalodisca)     
# A nice clean regular expression; \\d+ means at one or more digits. We're just stripping off any isoform designators.
Hvit.final.id <- gsub("_i", "", Hvit.final.id)
Hvit.final.id <- gsub(" ", "", Hvit.final.id)
Hvit.final.id <- data.frame(Hvit.final.id)
Hvit.final.id <- mutate(Hvit.final.id, gene.name = Hvit.gene.orthos$gene.name)
colnames(Hvit.final.id) <- c("Homalodisca", "gene.name")

rnaseqMatrix = HVdata[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=2,]

conditions = data.frame(conditions=factor(c(rep("HVAbd", 3), ### these are separate R scripts because there is missing data for H.vit.
                                            rep("HVEye", 2), 
                                            rep("HVLeg", 3), 
                                            rep("HVMeso", 3), 
                                            rep("HVOvi", 2), 
                                            rep("HVPro", 3),
                                            rep("HVWing2", 3),
                                            rep("HVWing3", 2))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)
HVdds = DESeq(ddsFullCountTable)

vst <- vst(HVdds)

## Reload data and start from here:
# load("Hvit_WingGenes_data.RData")

vst <- vst(HVdds)
vst <- assay(vst)
res = results(HVdds)
counts.dds <- data.frame(counts(HVdds, normalized=TRUE))

baseMeanAbd <- rowMeans(vst[,c(1,2,3)])
baseMeanEye <- rowMeans(vst[,c(4,5)])
baseMeanLeg <- rowMeans(vst[,c(6,7,8)])
baseMeanMeso <- rowMeans(vst[,c(9,10,11)])
baseMeanOvi <- rowMeans(vst[,c(12,13)])
baseMeanPro <- rowMeans(vst[,c(14,15,16)])
baseMeanWing2 <- rowMeans(vst[,c(17,18,19)])
baseMeanWing3 <- rowMeans(vst[,c(20,21)])

res = cbind(baseMeanAbd, baseMeanEye, baseMeanLeg, baseMeanMeso, baseMeanOvi, baseMeanPro, baseMeanWing2, baseMeanWing3, as.data.frame(res))
res = cbind(id=rownames(res), as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1

### write commands are commented out; uncomment if you want to write data to your machine.
#write.table(res, "Hvit_vst_means.tab",sep = "\t", quote = FALSE, row.names = res$id)
#save.image("Hvit_WingGenes_vstHeatmap.RData")
wingRes <- res %>% filter(id %in% Hvit.final.id$Homalodisca)
#write.table(wingRes, "Hvit_winggenes_vst_means.tab", sep="\t", quote = FALSE, row.names = wingRes$id)

#### heatmap.3.R file from github
source("https://raw.githubusercontent.com/trinityrnaseq/trinityrnaseq/devel/Analysis/DifferentialExpression/R/heatmap.3.R")

res2 <- wingRes[match(Hvit.final.id$Homalodisca, wingRes$id), ]
colnames(res2)[1] <- "Homalodisca"
gene.order <- inner_join(res2, Hvit.gene.orthos, by = "Homalodisca")
## Produced some doubles, probably from trying to deal with duplicated IDs in both species.
gene.order <- unique(gene.order) ## This keeps only each unique row, so if the Hvits don't match or the Ecars don't match it keeps it, which is what we want

## Now gene.order has the data in the right order
## So I can subset the data out of it
data <- as.matrix(gene.order[,c(8,9,7,2,5,4,6,3)])
gene.names <- gene.order$gene.name
rownames(data) <- gene.names
myheatcol = colorpanel(100, 'white','lavenderblush2','dodgerblue4') # sequential blues 
column_annotation <- c("violet", "purple", "cyan", "darkblue", "blue", "darkgreen", "lightgreen", "red") ## Body region color annotations
column_annotation <- as.matrix(t(column_annotation))

heatmap.3(data, dendrogram = "none", cexCol = .8, cexRow = .8, Rowv = FALSE, Colv = FALSE, 
          col = myheatcol, main="H.vitripennis VST(counts), wing candidate genes", 
          ColSideColors = column_annotation)


# pdf("HVWingGenes.vst-nonBlockwiseCounts-alphabeticalOrder2.pdf")
# heatmap.3(data, dendrogram = "none", cexCol = .8, cexRow = .8, Rowv = FALSE, Colv = FALSE, 
#           col = myheatcol, main="H.vitripennis VST(counts), wing candidate genes", 
#           ColSideColors = column_annotation)
# dev.off()


#~~~~~~~~ Uncomment these lines to save the data for reloading and rerunning ~~~~~~~~~~~#
# HVData <- data
# HVcolors <- myheatcol
# save(Hvit.final.id, HVdds, HVcolors, HVData, Hvit.gene.orthos, file="Hvit_WingGenes_data.RData")
