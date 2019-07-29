#!/bin/bash

source activate RNAseq_tetoolkit
snakemake -p --cores 48
conda deactivate
