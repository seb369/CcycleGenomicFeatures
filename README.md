# CcycleGenomicFeatures: Genomic features predict bacterial life history strategies in soil, as identified by metagenomic stable isotope probing

Author: Samuel Barnett

This repository contains the code used in the processing and analysis of the metagenomic DNA stable isotope probing study aiming to link genomic features to bacterial life history strategies in soil.

The paper to which this analysis belongs is titled "Genomic features predict bacterial life history strategies in soil, as identified by metagenomic stable isotope probing" (DOI:TBD) and can be cited as:

Barnett, S.E., Egan, R., Foster, B., Eloe-Fadrosh, E.A., and Buckley, D.H. (2022) Genomic features predict bacterial life history strategies in soil, as identified by metagenomic stable isotope probing. bioRxiv 507310.

## Directories

### Metagenome_processing
This directory contains code for some of the processing of the metagenomic data including:
* Aligning reads to contigs prior to binning
* Binning contigs into MAGs
* Annotating contigs using various tools (e.g., PROKKA, antismash)

### External_data_processing
This directory contains code for processing the data from external studies used to extend our findings beyon the scope of our metagenome. There where 7 studies used and this directory contains processing of whole genomes or MAGs from these studies, primarily annotation.

### Analysis
This directory contains the analysis of all data for this study. Analysis was performed in R using rMarkdown. Code can be found as rMarkdown files (.Rmd) or knitted markdown files (.md).

### publication_data_and_figures
Figures used in the publication made using the analysis described above.
