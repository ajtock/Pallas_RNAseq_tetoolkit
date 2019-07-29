#!/bin/bash

# Create conda environment for analysis of RNA-seq data
# using TEToolkit
conda env create --name RNAseq_tetoolkit --file environment.yaml

# DESeq2 was not included in original environment.yaml file,
# so add to RNAseq_tetoolkit environment
conda install --name RNAseq_tetoolkit --channel bioconda bioconductor-deseq2
