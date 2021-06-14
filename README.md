Curtiss Tfh RNA-seq
================
Chris Scharer

# Repository Info

Rscripts used for processing bulk RNA-seq data associated with the
following
publication:

Citation:

# Scripts for mapping, calculating gene coverage, and differential analysis

## mapping Folder

  - RNAseq.sample.manifest.txt: key file for all samples and contains
    metadata for each

  - RNAseq.pipeline.R: script for initial data organization, fastq
    qc/trimming, mapping, duplicate marking

## coverage Folder

  - RNAseq\_coverage\_Balbc.R: extracts coverage for all Entrez gene
    exons and normalizes data

## deseq2 folder

  - diff.glm.R: pairwise differential analysis with edgeR

## Common functions written for many data analysis routines

  - bisTools.R: common functions used for processing files
