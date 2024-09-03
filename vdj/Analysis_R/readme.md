This folder contains R scripts to analyse immune repertoire BCR/TCR data from cell Ranger. 

## Content: 
- alakazam.R: Analyse the repertoire diversity, gene usage and aminoacid properties. Input file: AIRR.tsv
- immunarch.R: Exploratory analysis of the repertoire. Gene usage, diversity, clonality, tracking clonotypes between samples ...
- integration_gex_bcr.R: Integrate BCR data to gene expression ddata. 1. Analysis of GEX data (seurat). 2. Integrate BCR.
- scRepertoire.R: Basics of the package + integration of GEX + B/TCR data.

STEPS: 
1. Analyse gene expression data with Seurat
2. Clonotype-only analysis (immune repertoire).
3. Integrate GEX and B/TCR. 
