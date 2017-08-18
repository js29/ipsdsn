# iPSC-derived sensory neurons: computational analyses
This repository contains code used for analysis of iPSC-derived sensory neurons, as described in this manuscript:

**Molecular and functional variation in iPSC-derived sensory neurons**
http://www.biorxiv.org/content/early/2017/01/06/095943

The code includes:
* Example pipeline for aligning RNA-seq reads, and quantifying expression
* Calling eQTLs and splice QTLs with FastQTL and with RASQUAL
* Producing plots included in the paper

## Directories
* **analysis** - scripts to analyse processed sensory neuron data
* **data** - input data for analyses (e.g. gene expression counts)
* **pipeline** - scripts to process raw data, e.g. alignment, quantification
* **utils** - utility scripts

Scripts within these directories may refer to files that are not part of the Github repository. Aligned BAM files are available from EGA (EGAD00001003145, 80 managed-access RNA-seq samples) or from ENA (accession ERP020576, 51 open-access RNA-seq samples). Similarly, ATAC-seq data for 8 and 23 managed and open-access samples are available. Summary statistics and gene expression counts are available at https://www.ebi.ac.uk/biostudies/studies/S-BSST16.

## Analysis
This directory includes mainly R scripts used for high-levels analyses, including producing figures in the paper.

## Data
This directory includes sample metadata.

## Pipelines
This directory includes scripts related aligning reads, computing gene expression counts, and calling expression QTLs and splice QTLs.

## Utils
This directory includes support scripts for pipelines, such as for submitting jobs to our computing cluster, running other software tools, and processing the results.
