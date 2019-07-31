#!/applications/R/R-3.5.0/bin/Rscript

# For each stranded paired-end library, reads were aligned to the
# TAIR10 reference genome using STAR version 2.7.0d.
# see /home/ajt200/analysis/20180928_Pallas_RNAseq_series1/fastq_pooled/snakemake_RNAseq_STAR/
# and /home/ajt200/analysis/20190215_Pallas_RNAseq_series2/fastq_pooled/snakemake_RNAseq_STAR/

# Transcript abundances were quantified using the TEcount script
# within TEToolkit version 2.0.3;
# see /home/ajt200/analysis/Pallas_RNAseq_tetoolkit/snakemake_RNAseq_tetoolkit

# Transcript-level estimates were summed by TEcount to derive a
# single expression estimate for each parent gene and TE identifier.

# This script:
# 1. Applies the regularized logarithm (rlog) transformation,
#    yielding approximately equal variances across mean expression estimates.
# 2. Calculates Euclidean distances between samples using the rlog-transformed data.
# 3. Generates principal component analysis and multi-dimensional scaling plots
#    for visualisation of sample-to-sample distances using the rlog-transformed data.

# R version 3.5.0
# DESeq2 version 1.22.2
# Note that this R version or later is required for DESeq2 version 1.16 or later,
# in which "the log2 fold change shrinkage is no longer default for the DESeq and nbinomWaldTest functions".
# DESeq2 version 1.16 introduces "a separate function lfcShrink, which performs log2 fold change shrinkage
# for visualization and ranking of genes." (see https://support.bioconductor.org/p/95695/ and
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#changes)


#Differentially expressed genes were identified using DESeq2 version 1.16.1 [74], using untransformed expression. Genes with more than one read across all samples within a contrast were retained. Additional filtering of genes with low mean read counts was automatically applied by DESeq2. For each contrast, differentially expressed genes with BH-adjusted P-values <0.01 were identified. Log2 fold change in gene expression was plotted against the mean of read counts normalized by library size for each gene in MA plots. A Bayesian method implemented in DESeq2 was used to moderate the log2 fold changes obtained for genes with low or variable expression levels. Up-regulated and down-regulated genes in taf4b-1 were evaluated for enrichment of genes up-regulated in wild type meiocytes compared to leaves (BH-adjusted P<0.01) using the hypergeometric distribution. Genes representing the intersection of those down-regulated, or up-regulated, in taf4b-1 (BH-adjusted P<0.01) and up-regulated in meiocytes (BH-adjusted P<0.01), were analyzed for gene ontology (GO) term enrichment. Gene sets were analyzed for over-representation of “biological process” GO terms relative to their representation among all genes in the TAIR10 annotation, using topGO (version 2.26.0) [75]. Significantly enriched terms were identified by applying the default topGO algorithm coupled with the Fisher’s exact test statistic (P≤0.05). 

# Usage:
# ./DESeq2_sample_exploration.R TEcount multi 'wt_RNAseq_Rep1,wt_RNAseq_Rep2,cmt3_RNAseq_Rep1,cmt3_RNAseq_Rep2,hta6_RNAseq_Rep1,hta6_RNAseq_Rep2,hta6_hta7_RNAseq_Rep1,hta6_hta7_RNAseq_Rep2,hta7_cmt3_RNAseq_Rep1,hta7_cmt3_RNAseq_Rep2,hta6_hta7_cmt3_RNAseq_Rep1,hta6_hta7_cmt3_RNAseq_Rep2' 'wt_Rep1,wt_Rep2,cmt3_Rep1,cmt3_Rep2,hta6_Rep1,hta6_Rep2,hta6hta7_Rep1,hta6hta7_Rep2,hta7cmt3_Rep1,hta7cmt3_Rep2,hta6hta7cmt3_Rep1,hta6hta7cmt3_Rep2'

#tool <- "TEcount"
#mode <- "multi"
#prefixes <- unlist(strsplit("wt_RNAseq_Rep1,wt_RNAseq_Rep2,cmt3_RNAseq_Rep1,cmt3_RNAseq_Rep2,hta6_RNAseq_Rep1,hta6_RNAseq_Rep2,hta6_hta7_RNAseq_Rep1,hta6_hta7_RNAseq_Rep2,hta7_cmt3_RNAseq_Rep1,hta7_cmt3_RNAseq_Rep2,hta6_hta7_cmt3_RNAseq_Rep1,hta6_hta7_cmt3_RNAseq_Rep2",
#                            split = ","))
#libNames <- unlist(strsplit("wt_Rep1,wt_Rep2,cmt3_Rep1,cmt3_Rep2,hta6_Rep1,hta6_Rep2,hta6hta7_Rep1,hta6hta7_Rep2,hta7cmt3_Rep1,hta7cmt3_Rep2,hta6hta7cmt3_Rep1,hta6hta7cmt3_Rep2",
#                            split = ","))

args <- commandArgs(trailingOnly = T)
tool <- args[1]
mode <- args[2]
prefixes <- unlist(strsplit(args[3],
                            split = ","))
libNames <- unlist(strsplit(args[4],
                            split = ","))

library(DESeq2)
print(packageVersion("DESeq2"))
#[1] ‘1.22.2’
library(dplyr)

outDir <- "DESeq2/"
plotDir <- paste0(outDir, "sample_exploration_plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load counts
counts <- lapply(seq_along(prefixes), function(x) {
  read.table(paste0(tool, "/", mode, "/",
                    tool, "_", mode, "_", prefixes[x], ".cntTable"),
                    header = T,
                    row.names = 1)
})

# Re-format so that featureID is column 1
counts <- lapply(seq_along(counts), function(x) {
  data_frame(featureID = rownames(counts[[x]]),
             counts = counts[[x]][,1])
})

# Merge columns
df <- counts %>% purrr::reduce(full_join,
                               by = "featureID")
df <- as.data.frame(df)
# Assign rownames using featureID column
rownames(df) <- df$featureID
# Remove featureID column
df <- df[,-1]
colnames(df) <- libNames

# Number of features before filtering
print("Features:")
print(nrow(df))
#[1] 37656

# Create table of sample names and conditions
sampleTable <- data.frame(sample = libNames,
                          condition = factor(sub(pattern = "_Rep\\d",
                                                 replacement = "",
                                                 x = libNames)))
rownames(sampleTable) <- colnames(df)
print(sampleTable)

# Retain only features that have more than a single read across all samples
df <- df[apply(X = df,
               MARGIN = 1,
               FUN = function(x) { sum(x) })
         > 1,]
print("Features with > 1 read count across all samples:")
print(nrow(df))
#[1] 37368

# Create DESeqDataSet (dds)
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = sampleTable,
                              design = ~ condition)

## The rlog and variance stabilizing transformations
# see http://www.bioconductor.org/help/workflows/rnaseqGene/#the-rlog-and-variance-stabilizing-transformations

rld <- rlog(dds, blind = FALSE)
print(head(assay(rld), 3))

vsd <- vst(dds, blind = FALSE)
print(head(assay(vsd), 3))

# Visualise the effect of transformation
library(dplyr)
library(ggplot2)
library(hexbin)

# For the log2 approach, estimate size factors to account for sequencing depth
# Sequencing-depth correction is done automatically for rlog and vst
dds <- estimateSizeFactors(dds)

trans_df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized = T)[,1:2]+1)) %>%
    mutate(transformation = "log2(normalized counts + 1)"),
  as_data_frame(assay(rld)[,1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[,1:2]) %>% mutate(transformation = "vst"))

#colnames(trans_df)[1:2] <- paste0(colnames(df)[1:2],
#                                  " transformed counts")

plot_transformed_counts <- ggplot(trans_df,
                                  aes(x = `wt_Rep1`,
                                      y = `wt_Rep2`)) +
                           geom_hex(bins = 80) +
                           coord_fixed() +
                           facet_grid(. ~ transformation) +
                           labs(fill = "Occurrences") +
                           theme_classic()
                           #theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
ggsave(plot_transformed_counts,
       file = paste0(plotDir,
                     colnames(df)[1], "_vs_",
                     colnames(df)[2], "_transformed_counts_log2countsPlus1_rlog_vst_genes_TEs.pdf"))


## Sample distances

# Sample distances using the rlog-transformed counts
sampleDists <- dist(t(assay(rld)))
print(sampleDists)

library(pheatmap)
library(RColorBrewer)

# Heatmap of sample distances using the rlog-transformed counts
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- as.vector(sampleTable$sample)
colnames(sampleDistMatrix) <- NULL
mycols <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf(paste0(plotDir, "sample_distances_heatmap_rlog_genes_TEs.pdf"),
    height = 5, width = 7.5, onefile = F)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = mycols)
dev.off()

# Sample distances using the Poisson Distance
library(PoiClaClu)
poisd <- PoissonDistance(t(counts(dds)))
print(poisd)
#Value of alpha used to transform data:  0.5353061
#This type of normalization was performed: mle
#Dissimilarity computed for  12  observations.

# Heatmap
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- as.vector(sampleTable$sample)
colnames(samplePoisDistMatrix) <- NULL
pdf(paste0(plotDir, "sample_distances_heatmap_Poisson_genes_TEs.pdf"), height = 5, width = 7.5, onefile = F)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = mycols)
dev.off()

# PCA plots for visualising sample-to-sample distances
PCAplot_rlog <- DESeq2::plotPCA(rld, intgroup = c("condition", "sample")) +
                  theme(plot.margin = grid::unit(c(0,0,0,0), "mm")) +
                  theme_classic() +
                  theme(panel.border = element_rect(colour = "black",
                                                    fill = NA,
                                                    size = 1)) +
                  coord_fixed()
ggsave(PCAplot_rlog,
       file = paste0(plotDir, "PCAplot_rlog_genes_TEs.pdf"),
       width = 20, height = 20, units = "cm")

PCAplot_rlog_data <- DESeq2::plotPCA(rld, intgroup = c("condition", "sample"),
                                     returnData = T)
print(PCAplot_rlog_data)
# Obtain percentage variance explained by PC1 and PC2 for plotting using ggplot2
percentVar <- round(100 * attr(PCAplot_rlog_data, "percentVar"))

PCAggplot_rlog <- ggplot(PCAplot_rlog_data,
                         aes(x = PC1, y = PC2,
                             shape = condition,
                             colour = sample)) +
                  geom_point(size = 3) +
                  scale_shape_manual(values = c(3, 8, 9, 15, 16, 17)) +
                  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
                  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
                  theme(plot.margin = grid::unit(c(0,0,0,0), "mm")) +
                  theme_classic() +
                  theme(panel.border = element_rect(colour = "black",
                                                    fill = NA,
                                                    size = 1)) +
                  coord_fixed()
ggsave(PCAggplot_rlog,
       file = paste0(plotDir, "PCAggplot_rlog_genes_TEs.pdf"),
       height = 20, width = 20, units = "cm")

# MDS (multi-dimensional scaling) plots for visualising sample-to-sample distances
# rlog-transformed counts
mds <- as.data.frame(colData(rld)) %>%
         cbind(cmdscale(sampleDistMatrix))
MDSplot_rlog <- ggplot(mds, aes(x = `1`, y = `2`,
                                shape = condition,
                                colour = sample)) +
                  geom_point(size = 3) +
                  scale_shape_manual(values = c(3, 8, 9, 15, 16, 17)) +
                  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
                  theme_classic() +
                  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
                  coord_fixed()
ggsave(MDSplot_rlog,
       file = paste0(plotDir, "MDSplot_rlog_genes_TEs.pdf"), height= 20, width = 20, units = "cm")

# Poisson distance
mdsPois <- as.data.frame(colData(dds)) %>%
             cbind(cmdscale(samplePoisDistMatrix))
MDSplot_Pois <- ggplot(mdsPois, aes(x = `1`, y = `2`,
                                    shape = condition,
                                    colour = sample)) +
                  geom_point(size = 3) +
                  scale_shape_manual(values = c(3, 8, 9, 15, 16, 17)) +
                  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
                  theme_classic() +
                  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
                  coord_fixed()
ggsave(MDSplot_Pois,
       file = paste0(plotDir, "MDSplot_Pois_genes_TEs.pdf"), height= 20, width = 20, units = "cm")

# feature clustering
# Clustering of 20 features with greatest rlog-transformed variance across samples
library(genefilter)
topVarfeatures <- head(order(rowVars(assay(rld)), decreasing = T), 20)
mat <- assay(rld)[topVarfeatures, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("sample", "condition")])
pdf(paste0(plotDir, "feature_clustering_rld_topVar20_genes_TEs.pdf"), height = 5, width = 7.5, onefile = F)
pheatmap(mat, annotation_col = anno)
dev.off()


################
#dds$groups = relevel(dds$groups,"CGroup")
#dds <- DESeq(dds)
#res <- results(dds,independentFiltering=F)
#write.table(res, file="TEtranscripts_out_gene_TE_analysis.txt", sep="\t",quote=F)
#resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) & (abs(res$log2FoldChange)> 0.000000)), ]
#write.table(resSig, file="TEtranscripts_out_sigdiff_gene_TE.txt",sep="\t", quote=F)
#
#df <- read.table("combined.cntTable",header=T,row.names=1)
## Assuming 2 treatment and 2 controls
#groups <- factor(c(rep("TGroup",2),rep("CGroup",2)))
#sampleInfo <- data.frame(groups,row.names=colnames(df))
#library(DESeq2, quietly=T)
#dds <- DESeqDataSetFromMatrix(countData = df, colData = sampleInfo, design = ~ groups)
#dds$condition = relevel(dds$groups,"CGroup")
#dds <- DESeq(dds)
#res <- results(dds,independentFiltering=F)
#write.table(res, file="pairedEnd_test_gene_TE_analysis.txt", sep="\t",quote=F)
#resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000), ]
#write.table(resSig, file="pairedEnd_test_sigdiff_gene_TE.txt",sep="\t", quote=F) 
