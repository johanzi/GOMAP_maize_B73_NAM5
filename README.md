# GOMAP_maize_B73_NAM5
GO term annotation of the maize B73 NAM5 assembly

# Objectives

I needed to perform GO enrichment analysis on genes targeted by sRNAs in maize before the time a B73 NAM5 GO term annotation was available on [MaizeGDB](https://download.maizegdb.org/GeneFunction_and_Expression/Pannzer_GO_Terms/). This annotation was performed using [PANNZER](http://ekhidna2.biocenter.helsinki.fi/sanspanz/), see reference paper [Törönen & Holm, 2022](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4193) GO term annotation available. I therefore decided to use the Gene Ontology Meta Annotator for Plants ([GOMAP]([https://bioinformapping.com/gomap/master/RUNNING.htm](https://bioinformapping.com/gomap/master/OVERVIEW.html)) pipeline.

Both PANNZER and GOMAP annotated 39,756 genes. However, GOMAP could generate 493,310 annotations vs 167,519 for PANNZER. Although GOMAP is recommended to be run on a HPC due to the sheer amount of calculation it needs, my HPC restriction did not allow me to run such a long job. It took my local workstation (10 cores of 2,2 GHz each) about one month to complete the annotation. I would therefore recommend anyone to save some time and energy here and use the provided annotation.

# Requisites

* UNIX-based computer with a decent amount of RAM and cores
* Singularity
* 


