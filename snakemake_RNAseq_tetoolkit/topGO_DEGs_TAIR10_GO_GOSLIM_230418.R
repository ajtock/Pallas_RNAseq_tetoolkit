#!/applications/R/R-3.5.0/bin/Rscript

#################################################################
# Analyse differentially expressed genes for GO term enrichment #
# relative to all TAIR10 genes on chromosomes                   #
#################################################################

# This script is based on a useful post by Avril Coghlan:
# http://avrilomics.blogspot.co.uk/2015/07/using-topgo-to-test-for-go-term.html

# Example usage:
# topGO_DEGs_TAIR10_GO_GOSLIM_230418.R hta7_v_wt FDR0.05_L2FC0.0 upReg genes BP 0.05 0.1

/home/ajt200/analysis/Pallas_RNAseq_tetoolkit/snakemake_RNAseq_tetoolkit/DESeq2/FDR0.05_L2FC0.0/hta7_v_wt/genes/res_hta7_v_wt_FDR0.05_L2FC0.0_chr_upRegSorted_genes_featureIDs.txt

TEcount multi 'cmt3_RNAseq_Rep1,cmt3_RNAseq_Rep2,wt_RNAseq_Rep1,wt_RNAseq_Rep2' 0.05 0.0 genes

#prefixes <- unlist(strsplit("cmt3_RNAseq_Rep1,cmt3_RNAseq_Rep2,wt_RNAseq_Rep1,wt_RNAseq_Rep2",
#                            split = ","))
#FDRnum <- 0.05
#FDRchar <- "0.05"
#L2FCnum <- 0.0
#L2FCchar <- "0.0"
#featureName <- "genes"

args <- commandArgs(trailingOnly = T)
mode <- args[2]
prefixes <- unlist(strsplit(args[3],
                            split = ","))
FDRnum <- as.numeric(args[4])
FDRchar <- as.character(args[4])
L2FCnum <- as.numeric(args[5])
L2FCchar <- as.character(args[5])
featureName <- args[6]

contrast <- paste0(geno1, "_v_", geno2)

#options(echo=TRUE) # if you want to see commands in output file
args <- commandArgs(trailingOnly = TRUE)
ont <- args[1]
target <- args[2]
sigLevel <- args[3]
FDRchar <- as.character(args[4])

inDir <- paste0("./FDR", FDRchar, "/")

# load topGO
suppressMessages(library(topGO))
suppressMessages(library(Rgraphviz))

genesetGO <- function(target) {
  # Read in GO annotations for TAIR10 genes to define "gene universe"
  geneID2GO <- readMappings(file = paste0("/projects/ajt200/TAIR10/TAIR10_ATH_GO_GOSLIM_030418_geneID_GOann_reshaped_noC_noM.txt"))
  geneUniverse <- names(geneID2GO)

  # Define list of genes of interest; file should contain a single column of gene identifiers
  genesOfInterest <- as.character(read.table(paste0(inDir, target))$x)

  # Specify where genes of interest appear in the gene universe vector
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse

  # Build GOdata object in topGO
  capture.output(GOdata <- new("topGOdata", description = "Differentially expressed genes", ontology = ont,
                               allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO),
                 file="/dev/null")

  # Access list of genes of interest
  #sg <- sigGenes(GOdata)
  #print(str(sg))
  #print(numSigGenes(GOdata))

  # Run Fisher's exact tests to determine GO term enrichment
  capture.output(resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher"),
                 file="/dev/null")
  capture.output(resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher"),
                 file="/dev/null")
  capture.output(resultTopGO <- runTest(GOdata, algorithm = "weight01", statistic = "fisher"),
                 file="/dev/null")

  # Count number of results where weight01 gives a P-value <= sigLevel (args[3])
  mySummary <- summary(attributes(resultTopGO)$score <= as.numeric(sigLevel))
  numSignif <- as.integer(mySummary[[3]])

  # List significant results and write to file
  capture.output(enrichRes <- GenTable(GOdata, classicFisher = resultClassic,
                                       elimFisher = resultElim,
                                       topGOFisher = resultTopGO,
                                       orderBy = "topGOFisher",
                                       ranksOf = "elimFisher", topNodes = numSignif),
                 file="/dev/null")

  # WARNING: ASSUMES INPUT FILE HAS A 3-LETTER EXTENSION
  basename <- basename(target)
  len <- nchar(basename)
  basename <- substr(basename, 1, len-4)
  
  out_name <- paste(basename, "GO", ont, "enrichment.tsv", sep="_")
  folder <- paste0(inDir, "GO")
  system(paste0("[ -d ", folder, " ] || mkdir ", folder))
  folder2 <- paste0(folder, "/", basename, "_GO_", ont)
  system(paste0("[ -d ", folder2, " ] || mkdir ", folder2))

  capture.output(write.table(enrichRes, file = file.path(folder, out_name), sep = "\t",
                             row.names = FALSE, col.names = TRUE, quote = FALSE),
                 file="/dev/null")
  
  # Visualise the positions of the top 5 statistically significant GO terms in the GO hierarchy
  out_name2 <- paste(basename, "GO", ont, "enrichment", sep="_")
  printGraph(GOdata, resultTopGO, firstSigNodes = 5,
             fn.prefix = file.path(folder, out_name2), useInfo = "all", pdfSW = TRUE)

  # Extract gene IDs annotated with significantly enriched GO terms
  myTerms <- enrichRes$GO.ID
  myGenes <- genesInTerm(GOdata, myTerms)
  for(i in 1:length(myTerms)) {
    myTerm <- myTerms[i]
    myGenesForTerm <- myGenes[myTerm][[1]]
    myFactor <- myGenesForTerm %in% genesOfInterest
    myGenesForTermT <- myGenesForTerm[myFactor == TRUE]
    myGenesForTermT <- paste(myGenesForTermT, collapse = ",")
    myGenesForTermT <- paste(myTerm, myGenesForTermT, sep = "\t")
    out_name3 <- paste0(basename, "_GO_", ont, "_enrichment_", myTerm, ".txt")  
    write(myGenesForTermT, file = file.path(folder2, out_name3))
  }
}

# Apply degGO() function to target (differentially expressed genes file)
genesetGO(target = target)


sessionInfo()

