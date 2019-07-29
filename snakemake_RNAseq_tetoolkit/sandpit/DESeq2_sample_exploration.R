#!/applications/R/R-3.5.0/bin/Rscript

# For each stranded paired-end library, reads were aligned to the
# TAIR10 reference genome using STAR version 2.7.0d.
# see /home/ajt200/analysis/20180928_Pallas_RNAseq_series1/fastq_pooled/snakemake_RNAseq_STAR/
# and /home/ajt200/analysis/20190215_Pallas_RNAseq_series2/fastq_pooled/snakemake_RNAseq_STAR/

# Transcript abundances were quantified using the TEcounts script
# within TEToolkit version 2.0.3;
# see /home/ajt200/analysis/Pallas_RNAseq_tetoolkit/snakemake_RNAseq_tetoolkit

# Transcript-level estimates were summed by TEcounts to derive a
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

library(DESeq2)
print(packageVersion("DESeq2"))
#[1] ‘1.22.2’

# Load count data
data <- read.table("TEtranscripts_out.cntTable",
                   header = TRUE,
                   row.names = 1)
print("Features:")
print(nrow(data))
print(colnames(data))
colnames(data) <- c(
                    "hta6 Rep1",
                    "hta6 Rep2",
                    "wt Rep1",
                    "wt Rep2"
                   )
sampleTable <- data.frame(sample = colnames(data),
                          condition = factor(rep.int(c(
                                                       "hta6",
                                                       "wt"
                                                      ),
                                                     times = c(2, 2))))
print(sampleTable)


# Retain only features that have more than a single read across all samples
data <- data[apply(X = data,
                   MARGIN = 1,
                   FUN = function(x) { sum(x) })
             > 1,]
print("Features with > 1 read count across all samples:")
print(nrow(data))

# Create DESeqDataSet (dds)
dds <- DESeqDataSetFromMatrix(countData = data,
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

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized = T)[,1:2]+1)) %>%
    mutate(transformation = "log2(normalized counts + 1)"),
  as_data_frame(assay(rld)[,1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[,1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("Sample 1 transformed counts",
                       "Sample 2 transformed counts")

plot_transformed_counts <- ggplot(df,
                                  aes(x = `Sample 1 transformed counts`,
                                      y = `Sample 2 transformed counts`)) +
                           geom_hex(bins = 80) +
                           coord_fixed() +
                           facet_grid(. ~ transformation) +
                           labs(fill = "Occurrences") +
                           theme_classic()
                           #theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
ggsave(plot_transformed_counts,
       file = paste0(plotDir,
                     "Sample1_vs_Sample2_transformed_counts_log2countsPlus1_rlog_vst_genes.pdf"))

#groups <- factor(c(rep("hta6", 2),
#                   rep("wt", 2)))
#sampleInfo <- data.frame(groups,
#                         row.names = colnames(data))
#dds2 <- DESeqDataSetFromMatrix(countData = data,
#                               colData = sampleInfo,
#                               design = ~ groups)

dds$groups = relevel(dds$groups,"CGroup")
dds <- DESeq(dds)
res <- results(dds,independentFiltering=F)
write.table(res, file="TEtranscripts_out_gene_TE_analysis.txt", sep="\t",quote=F)
resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) & (abs(res$log2FoldChange)> 0.000000)), ]
write.table(resSig, file="TEtranscripts_out_sigdiff_gene_TE.txt",sep="\t", quote=F)

data <- read.table("combined.cntTable",header=T,row.names=1)
# Assuming 2 treatment and 2 controls
groups <- factor(c(rep("TGroup",2),rep("CGroup",2)))
sampleInfo <- data.frame(groups,row.names=colnames(data))
library(DESeq2, quietly=T)
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds$condition = relevel(dds$groups,"CGroup")
dds <- DESeq(dds)
res <- results(dds,independentFiltering=F)
write.table(res, file="pairedEnd_test_gene_TE_analysis.txt", sep="\t",quote=F)
resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000), ]
write.table(resSig, file="pairedEnd_test_sigdiff_gene_TE.txt",sep="\t", quote=F) 
