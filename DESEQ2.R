library(DESeq2)
library("pheatmap", lib="/scratch/ad163/Manuka2/Project_Aissatou/Count_htseq/ALL_MANUKA/pheatmap/")

options(bitmapType="cairo")

directory<-"/scratch/xx12/CBFa/Results_February_2015/P_Percipalle_14_02/Analysis/Counts/"

sampleFiles<-grep(".txt",list.files(directory),value=TRUE)
sampleFiles

sampleN <- sub("(*)_*.txt","\\1",sampleFiles)

sampleCondition <- factor(c(rep("het", 3),rep("kno", 4),rep("wt", 4)))
sampleCondition

sampleTable <- data.frame(sampleName = sampleN, fileName = sampleFiles, condition = sampleCondition)
sampleTable


dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = directory,design= ~ condition,ignoreRank=FALSE)
dds

dds <- DESeq(dds)

rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)

write.table(rlogMat, "/scratch/xx12/CBFa/Results_February_2015/P_Percipalle_14_02/Analysis/Counts/DESeq2/deseq2_all_samples_rlog.csv", sep="\t")
write.table(vstMat, "/scratch/xx12/CBFa/Results_February_2015/P_Percipalle_14_02/Analysis/Counts/DESeq2/deseq2_all_samples_vst.csv", sep="\t")


sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png("/scratch/xx12/CBFa/Results_February_2015/P_Percipalle_14_02/Analysis/Counts/DESeq2/deseq2_all_samples_rlog_heatmap.png")
rlogHeat <- pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)
rlogHeat
dev.off()


sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png("/scratch/xx12/CBFa/Results_February_2015/P_Percipalle_14_02/Analysis/Counts/DESeq2/deseq2_all_samples_vst_heatmap.png")
vstHeat <- pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)
vstHeat
dev.off()



ChetVSCkno <- results(dds,contrast=c("condition","het","kno"))
write.table(ChetVSCkno, "/scratch/xx12/CBFa/Results_February_2015/P_Percipalle_14_02/Analysis/Counts/DESeq2/hetVSkno.csv", sep="\t")
png("/scratch/xx12/CBFa/Results_February_2015/P_Percipalle_14_02/Analysis/Counts/DESeq2/hetVSkno_MAplot.png")
plotMA(ChetVSCkno, main="het_VSC_kno", ylim=c(-8,8))
dev.off()

ChetVSCwt <- results(dds,contrast=c("condition","het","wt"))
write.table(ChetVSCwt, "/scratch/xx12/CBFa/Results_February_2015/P_Percipalle_14_02/Analysis/Counts/DESeq2/hetVSwt.csv", sep="\t")
png("/scratch/xx12/CBFa/Results_February_2015/P_Percipalle_14_02/Analysis/Counts/DESeq2/hetVSwt_MAplot.png")
plotMA(ChetVSCwt, main="het_VSC_wt", ylim=c(-8,8))
dev.off()

CknoVSCwt <- results(dds,contrast=c("condition","kno","wt"))
write.table(CknoVSCwt, "/scratch/xx12/CBFa/Results_February_2015/P_Percipalle_14_02/Analysis/Counts/DESeq2/knoVSwt.csv", sep="\t")
png("/scratch/xx12/CBFa/Results_February_2015/P_Percipalle_14_02/Analysis/Counts/DESeq2/knoVSwt_MAplot.png")
plotMA(CknoVSCwt, main="kno_VSC_wt", ylim=c(-8,8))
dev.off()

