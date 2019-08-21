#!/applications/R/R-3.3.2/bin/Rscript

########################################
# Analyse genes for GO term enrichment #
# relative to all TAIR10 genes         #
########################################

# This script is based on a useful post by Avril Coghlan:
# http://avrilomics.blogspot.co.uk/2015/07/using-topgo-to-test-for-go-term.html

# Doesn't work with mclapply or dopar

# Example usage:
# /applications/R/R-3.3.2/bin/Rscript ./topGO_DEGs_TAIR10_GO_GOSLIM_030418_incl_weight.R BP 0.05 up 'FDR0.05_L2FC0.0'

#ont <- "BP"
#sigLevel <- 0.05
#reg <- "up"
#DEsigLevel <- "FDR0.05_L2FC0.0"

#options(echo=TRUE) # if you want to see commands in output file
args <- commandArgs(trailingOnly = TRUE)
ont <- args[1]
sigLevel <- args[2]
reg <- args[3]
DEsigLevel <- args[4]

suppressMessages(library(topGO))
suppressMessages(library(Rgraphviz))

contrasts <- c(
               "cmt3_v_hta6_hta7",
               "cmt3_v_hta6_hta7_cmt3",
               "cmt3_v_wt",
               "hta6_hta7_cmt3_v_wt",
               "hta6_hta7_v_wt",
               "hta6_v_hta6_hta7_cmt3",
               "hta6_v_hta7_cmt3",
               "hta6_v_wt",
               "hta7_cmt3_v_wt",
               "hta7_v_cmt3",
               "hta7_v_hta6_hta7_cmt3",
               "hta7_v_hta7_cmt3",
               "hta7_v_wt"
              )
targets <- paste0(DEsigLevel, "/", contrasts, "/",
                  "genes/",
                  "res_", contrasts, "_",
                  DEsigLevel, "_chr_",
                  reg, "RegSorted_genes_featureIDs.txt")

# Read in GO annotations for TAIR10 genes to define "gene universe"
geneID2GO <- readMappings(file = paste0("/projects/ajt200/TAIR10/TAIR10_ATH_GO_GOSLIM_030418_geneID_GOann_reshaped_noC_noM.txt"))
geneUniverse <- names(geneID2GO)

genesetGO <- function(target) {
  # Define list of genes of interest; file should contain a single column of gene identifiers
  genesOfInterest <- as.character(read.table(target)$V1)
  
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
  capture.output(resultWeight <- runTest(GOdata, algorithm = "weight", statistic = "fisher"),
                 file="/dev/null")
  capture.output(resultTopGO <- runTest(GOdata, algorithm = "weight01", statistic = "fisher"),
                 file="/dev/null")
  
  # Count number of results where weight01 gives a P-value <= sigLevel (args[3])
  mySummary <- summary(attributes(resultTopGO)$score <= as.numeric(sigLevel))
  numSignif <- as.integer(mySummary[[3]])
  
  # List significant results and write to file
  capture.output(enrichRes <- GenTable(GOdata,
                                       classicFisher = resultClassic,
                                       elimFisher = resultElim,
                                       weightFisher = resultWeight,
                                       topGOFisher = resultTopGO,
                                       orderBy = "topGOFisher",
                                       ranksOf = "elimFisher",
                                       topNodes = numSignif),
                 file="/dev/null")
  
  # WARNING: ASSUMES INPUT FILE HAS A 3-LETTER EXTENSION
  basename <- basename(target)
  len <- nchar(basename)
  basename <- substr(basename, 1, len-4)
  
  out_name <- paste(basename, "GO", ont, "enrichment.tsv", sep="_")
  folder <- paste0(dirname(target), "/GO")
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

# Apply genesetGO() function to each target (differentially expressed genes file)
lapply(seq_along(targets), function(x) {
  genesetGO(target = targets[x])
})
