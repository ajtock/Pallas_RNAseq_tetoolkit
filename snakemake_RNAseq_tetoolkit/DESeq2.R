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

# This script performs differential expression analysis using
# DESeq2 version 1.22.2

# R version 3.5.0
# DESeq2 version 1.22.2
# Note that this R version or later is required for DESeq2 version 1.16 or later,
# in which "the log2 fold change shrinkage is no longer default for the DESeq and nbinomWaldTest functions".
# DESeq2 version 1.16 introduces "a separate function lfcShrink, which performs log2 fold change shrinkage
# for visualization and ranking of genes." (see https://support.bioconductor.org/p/95695/ and
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#changes)

# Differentially expressed genes were identified using DESeq2 version 1.22.2, using untransformed expression value in accordance with DESeq2 model-fitting assumptions. Genes with more than one read across all samples within a contrast were retained. Additional filtering of genes with low mean read counts was automatically applied by DESeq2. For each contrast, differentially expressed genes with BH-adjusted P-values <0.05 were identified. Log2 fold change in gene expression was plotted against the mean of read counts normalized by library size for each gene in MA plots. A Bayesian method implemented in DESeq2 was used to moderate the log2 fold changes obtained for genes with low or variable expression levels. Up-regulated and down-regulated genes in each mutant were evaluated for enrichment of genes up-regulated in wild type meiocytes compared to leaves (BH-adjusted P<0.01) using the hypergeometric distribution. Genes representing the intersection of those down-regulated, or up-regulated, in each mutant (BH-adjusted P<0.05) and up-regulated in meiocytes (BH-adjusted P<0.01), were analyzed for gene ontology (GO) term enrichment. Gene sets were analyzed for over-representation of “biological process” GO terms relative to their representation among all genes in the TAIR10 annotation, using topGO (version 2.26.0) [75]. Significantly enriched terms were identified by applying the default topGO algorithm coupled with the Fisher’s exact test statistic (P≤0.05). 

# Usage:
# ./DESeq2.R TEcount multi 'cmt3_RNAseq_Rep1,cmt3_RNAseq_Rep2,wt_RNAseq_Rep1,wt_RNAseq_Rep2' 'cmt3_Rep1,cmt3_Rep2,wt_Rep1,wt_Rep2' 0.05 0.0 genes

#tool <- "TEcount"
#mode <- "multi"
#prefixes <- unlist(strsplit("cmt3_RNAseq_Rep1,cmt3_RNAseq_Rep2,wt_RNAseq_Rep1,wt_RNAseq_Rep2",
#                            split = ","))
#libNames <- unlist(strsplit("cmt3_Rep1,cmt3_Rep2,wt_Rep1,wt_Rep2",
#                            split = ","))
#FDRnum <- 0.05
#FDRchar <- "0.05"
#L2FCnum <- 0.0
#L2FCchar <- "0.0"
#featureName <- "genes"

args <- commandArgs(trailingOnly = T)
tool <- args[1]
mode <- args[2]
prefixes <- unlist(strsplit(args[3],
                            split = ","))
libNames <- unlist(strsplit(args[4],
                            split = ","))
FDRnum <- as.numeric(args[5])
FDRchar <- as.character(args[5])
L2FCnum <- as.numeric(args[6])
L2FCchar <- as.character(args[6])
featureName <- args[7]

geno1 <- sub(pattern = "_Rep\\d",
             replacement = "",
             x = libNames[1])
geno2 <- sub(pattern = "_Rep\\d",
             replacement = "",
             x = libNames[3])
contrast <- paste0(geno1, "_v_", geno2)

library(DESeq2)
print(packageVersion("DESeq2"))
#[1] ‘1.22.2’
library(dplyr)

featureDir <- paste0("DESeq2/", featureName, "/") 
sigDir <- paste0(featureDir, "FDR", FDRchar, "_L2FC", L2FCchar, "/")
outDir <- paste0(sigDir, contrast, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", featureDir, " ] || mkdir ", featureDir))
system(paste0("[ -d ", sigDir, " ] || mkdir ", sigDir))
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

# Retain rows corresponding to featureName
if(featureName == "genes") {
  df <- df[grepl(pattern = "^AT\\wG\\d+",
                 x = df$featureID),]
} else if(featureName == "TEs") {
  df <- df[!(grepl(pattern = "^AT\\wG\\d+",
                   x = df$featureID)),]
} else {
  stop("featureName is neither genes nor TEs")
}

# Assign rownames using featureID column
rownames(df) <- df$featureID
# Remove featureID column
df <- df[,-1]
colnames(df) <- libNames

# Number of features before filtering
print("Features:")
print(nrow(df))
#[1] 37336

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
#[1] 36344

# Create DESeqDataSet (dds)
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = sampleTable,
                              design = ~ condition)

# Fit model
dds_fit <- DESeq(dds)
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
res <- results(object = dds_fit,
               contrast = c("condition",
                            geno1,
                            geno2),
               alpha = FDRnum,
               lfcThreshold = L2FCnum,
               independentFiltering = TRUE,
               pAdjustMethod = "BH")
print(mcols(res, use.names = T))
print(summary(res))
print(paste0("DE ", featureName, ":"))
print(table((res$padj < FDRnum &
             res$log2FoldChange > L2FCnum) |
            (res$padj < FDRnum &
             res$log2FoldChange < (L2FCnum)*-1)))
print(paste0(featureName, " with non-NA FDR values"))
print(sum(!is.na(res$padj)))


## MA-plots
# These show log2(fold change) in gene expression
# for a variable (e.g., treated vs untreated)
# against the mean of normalised counts for the samples

# DESeq2 provides for use of a "Bayesian procedure to moderate
# (or "shrink") log2 fold changes from genes with very low counts
# and highly variable counts"
# Before generating an MA-plot, the lfcShrink() function should
# be used to moderate the log2(fold changes) for the contrast
res_lfcShrink <- lfcShrink(dds_fit,
                           contrast = c("condition",
                                        geno1,
                                        geno2),
                           res = res)
pdf(paste0(plotDir, "MAplot_", contrast,
           "_FDR", FDRchar, "_L2FC", L2FCchar,
           "_lfcShrink_", featureName, ".pdf"),
    height= 5, width = 6)
par(mfrow = c(1, 1))
par(mar = c(4.1, 4.1, 2.1, 2.1))
plotMA(res_lfcShrink,
       ylim = c(min(res_lfcShrink$log2FoldChange),
                max(res_lfcShrink$log2FoldChange)),
       xlab = "Mean normalized expression",
       ylab = bquote("Log"[2]*"("*.(geno1)*"/"*.(geno2)*")"),
       colSig = "darkorange1",
       colNonSig = "dimgrey",
       main = bquote(.(geno1)*" vs "*.(geno2)*
                     " ("*.(featureName)*
                     " with FDR < "*.(FDRchar)*
                     " and L2FC > "*.(L2FCchar)*" )"))
if(L2FCnum > 0) {
  abline(h = c(L2FCnum, (L2FCnum)*-1),
         lwd = 1.5, lty = 2, col = "red")
}
dev.off()


## Exporting results
# Exclude mitochondrial and chloroplast features
res_chr <- res[!grepl(pattern = "^ATM", x = rownames(res)) &
               !grepl(pattern = "^ATC", x = rownames(res)),]
# Subset results table to significant (padj < FDRnum)
# down-regulated genes (log2FoldChange < (L2FC)*-1)
res_chr_downReg <- res_chr[!is.na(res_chr$padj) &
                           res_chr$padj < FDRnum &
                           res_chr$log2FoldChange < (L2FCnum)*-1,]
# OR
#res_chr_downReg2 <- subset(res_chr,
#                           padj < FDRnum &
#                           log2FoldChange < (L2FCnum)*-1)

# Sort by increasing log2 fold change estimate
# (genes with strongest down-regulation at top)
res_chr_downRegSorted <- res_chr_downReg[order(res_chr_downReg$log2FoldChange),]
res_chr_downRegSorted_featureIDs <- rownames(res_chr_downRegSorted)
write.table(res_chr_downRegSorted,
            file = paste0(outDir, "res_", contrast,
                          "_FDR", FDRchar, "_L2FC", L2FCchar,
                          "_chr_downRegSorted_", featureName, ".tsv"),
            sep = "\t", quote = F)
write.table(res_chr_downRegSorted_featureIDs,
            file = paste0(outDir, "res_", contrast,
                          "_FDR", FDRchar, "_L2FC", L2FCchar,
                          "_chr_downRegSorted_", featureName, "_featureIDs.txt"),
            quote = F, row.names = F, col.names = F)

# Subset results table to significant (padj < FDRnum)
# up-regulated genes (log2FoldChange > L2FC)
res_chr_upReg <- res_chr[!is.na(res_chr$padj) &
                         res_chr$padj < FDRnum &
                         res_chr$log2FoldChange > L2FCnum,]
# OR
#res_chr_upReg2 <- subset(res_chr,
#                         padj < FDRnum &
#                         log2FoldChange > L2FCnum)

# Sort by decreasing log2 fold change estimate
# (genes with strongest up-regulation at top)
res_chr_upRegSorted <- res_chr_upReg[order(res_chr_upReg$log2FoldChange,
                                           decreasing = T),]
res_chr_upRegSorted_featureIDs <- rownames(res_chr_upRegSorted)
write.table(res_chr_upRegSorted,
            file = paste0(outDir, "res_", contrast,
                          "_FDR", FDRchar, "_L2FC", L2FCchar,
                          "_chr_upRegSorted_", featureName, ".tsv"),
            sep = "\t", quote = F)
write.table(res_chr_upRegSorted_featureIDs,
            file = paste0(outDir, "res_", contrast,
                          "_FDR", FDRchar, "_L2FC", L2FCchar,
                          "_chr_upRegSorted_", featureName, "_featureIDs.txt"),
            quote = F, row.names = F, col.names = F)

