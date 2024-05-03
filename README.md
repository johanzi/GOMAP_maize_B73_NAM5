# GOMAP_maize_B73_NAM5
GO term annotation of the maize B73 NAM5 assembly

# Objectives

I needed to perform GO enrichment analysis on genes targeted by sRNAs in maize before the time a B73 NAM5 GO term annotation was available on [MaizeGDB](https://download.maizegdb.org/GeneFunction_and_Expression/Pannzer_GO_Terms/) in 2022. This annotation was performed using [PANNZER](http://ekhidna2.biocenter.helsinki.fi/sanspanz/), see reference paper [Törönen & Holm, 2022](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4193). I therefore decided to use the Gene Ontology Meta Annotator for Plants ([GOMAP]([https://bioinformapping.com/gomap/master/RUNNING.htm](https://bioinformapping.com/gomap/master/OVERVIEW.html)) pipeline.

Both PANNZER and GOMAP annotated 39,756 genes. However, GOMAP could generate 493,310 annotations vs 167,519 for PANNZER. Although GOMAP manual recommends running the pipeline on a high-performance cluster (HPC) due to the sheer amount of calculation it needs, the setting of my HPC provider did not allow me to run such a long job. It took my local workstation (10 quad-core processors of 2,2 GHz each) about one month to complete the annotation. To save time and energy for others, I provide here the GOMAP annotation for B73 NAM5 with its step-by-step.

# Requisites

* UNIX-based computer with a decent amount of RAM and cores (mine had 10 quad-core processors (40 cores) of 2.2 GHz each, 63 Gb of RAM, 4 Gb GPU, OS Linux Mint 19.1 Cinnamon)
* Singularity (version 3.6.3 used here
* GOMAP (version 1.3.5 used here)

# Install Singularity and go

I found this step to be the most tedious somehow. I had to use several come-arounds to install the software on my workstation. I recommend you first try the installation method from [syslabs](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) and see if you need to go through the hard way (below) or not. I put the script I used here but you may try an easier way and you will have to adjust some variables.

```
sudo apt-get update && \
sudo apt-get install -y build-essential \
libseccomp-dev pkg-config squashfs-tools cryptsetup

sudo rm -r /usr/local/go

export VERSION=1.13.15 OS=linux ARCH=amd64  # change this as you need

wget -O /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz https://dl.google.com/go/go${VERSION}.${OS}-${ARCH}.tar.gz && \
sudo tar -C /usr/local -xzf /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz

echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
source ~/.bashrc

curl -sfL https://install.goreleaser.com/github.com/golangci/golangci-lint.sh |
sh -s -- -b $(go env GOPATH)/bin v1.21.0

mkdir -p ${GOPATH}/src/github.com/sylabs && \
cd ${GOPATH}/src/github.com/sylabs && \
git clone https://github.com/sylabs/singularity.git && \
cd singularity

git checkout v3.6.3

cd ${GOPATH}/src/github.com/sylabs/singularity && \
./mconfig && \
cd ./builddir && \
make && \
sudo make install

singularity version
```

# Install GOMAP

Once Singularity and go are installed properly, GOMAP should be a breeze to install.

```
# Clone GOMAP git repository
git clone https://github.com/Dill-PICL/GOMAP-singularity.git 
git checkout v1.3.5

cd GOMAP-singularity/

# Add this to your ~/.bashrc or run the line in the terminal
export GOMAP_LOC="/path/to/GOMAP-singularity/"

# Check if GOMAP runs
./test.sh
```


# Protein sequences of B73 NAM5

Download fasta file containing all annotated peptides from Maize GDB.

```
# Download fasta file
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa

# Check number of entries
grep ">" Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa | wc -l
72539

# Remove all asterisks ending sequences (would create a bug in GOMAP)
sed -i 's/\*//' Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa
```

# Create config file for B73 NAM5

GOMAP requires a configuration file called `min-config.yml` to run. I downloaded a template from the gomap developers and tweaked it a bit.
```
# Get config file for that job
wget --no-check-certificate https://bioinformapping.com/gomap/master/_static/min-config.yml min-config.yml.1
```

I changed the variables to get this:

```
input:
  #input fasta file name
  fasta: Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa
  # output file basename
  basename: maize_B73_NAM5
  #input NCBI taxonomy id
  taxon: "4577"
  # Name of the species
  species: "Zea mays"
  # Email is mandatory
  email: your_email@gmail.com
  #Number of CPUs used for tools
  cpus: 10
  #Whether mpi should be used (mpich-3.2.1  is default)
  mpi: False
  #what the name of the 
  tmpdir: "/tmp"
```

Note some important things: 
* Location of the fasta file, which should be in the same directory as the min-config.yml
* basename, this will be the prefix (maize_B73_NAM5) of the output gaf (GO annotation file format) file, here `maize_B73_NAM5.aggregate.gaf`
* Your email, needed to run one of the tools (Argot GO)
* cpus (=cores), depends on your workstation, the more, the faster your job will be. Plan anyway a good month to let your machine run if you don't have access to an HPC

# GOMAP run

GOMAP runs four independent steps:

1. seqsim
2. domain
3. fanngo
4. mixmeth-blast

Run the 4 first steps in parallel, each using 10 cores (40 cores in total available on my workstation)

After 1-4 steps are finished, the last three steps need to be run consecutively
5. mixmeth-preproc
6. mixmeth
7. aggregate

Here are the commands:

```
cd /path/to/GOMAP-singularity/

./run-GOMAP-SINGLE.sh --step=seqsim --config=B73_NAM5/min-config.yml

./run-GOMAP-SINGLE.sh --step=domain --config=B73_NAM5/min-config.yml

./run-GOMAP-SINGLE.sh --step=fanngo --config=B73_NAM5/min-config.yml

./run-GOMAP-SINGLE.sh --step=mixmeth-blast --config=B73_NAM5/min-config.yml

./run-GOMAP-SINGLE.sh --step=mixmeth-preproc --config=B73_NAM5/min-config.yml

./run-GOMAP-SINGLE.sh --step=mixmeth --config=B73_NAM5/min-config.yml

./run-GOMAP-SINGLE.sh --step=aggregate --config=B73_NAM5/min-config.yml
```

# Post-processing

The output file will look like this:

```
head maize_B73_NAM5.aggregate.gaf
!gaf-version:2.0
!db     db_object_id    db_object_symbol        qualifier       term_accession  db_reference    evidence_code   with    aspect  db_object_name  db_object_synonym       db_object_type  taxon   date    assigned_by     annotation_extensiongene_product_form_id
GOMAP   Zm00001eb000010_P001    Zm00001eb000010_P001            GO:0003690      GOMAP:0000      IEA     SMART:SM00733   F                       gene    taxon:4577      20221204        GOMAP-v1.3.9
GOMAP   Zm00001eb000010_P001    Zm00001eb000010_P001            GO:0003727      GOMAP:0000      IEA     RBH:I1IST9      F                       gene    taxon:4577      20221129        GOMAP-v1.3.9
```

Note that each protein isoforms are annotated but most people usually consider only genes when performing GO term enrichment analyses. Therefore, one must collapse this annotation at the gene level:

```
# Get protein isoform and GO term
# Remove 2 first lines (headers)
sed '1,2d' maize_B73_NAM5.aggregate.gaf | cut -f2 - | cut -d_ -f1 - > maize_B73_NAM5.aggregate.gaf.genes
sed '1,2d' maize_B73_NAM5.aggregate.gaf | cut -f5 - > maize_B73_NAM5.aggregate.gaf.go

paste maize_B73_NAM5.aggregate.gaf.go maize_B73_NAM5.aggregate.gaf.genes | sort | uniq > maize_B73_NAM5.aggregate.GO_gene.gaf

# Compare number of lines once protein isoforms are merged
wc -l maize_B73_NAM5.aggregate.gaf
846030 maize_B73_NAM5.aggregate.gaf

wc -l maize_B73_NAM5.aggregate.GO_gene.gaf
493310 maize_B73_NAM5.aggregate.GO_gene.gaf
```

Note the difference, with 846,030 annotations when considering protein isoforms, and 493,310 when collapsed at the gene level. If different protein isoforms of a same gene have different GO terms, the information is kept at the gene level so this way is OK for most people I think.

# GO term enrichment analysis

Now, one can use this gaf file in different softwares performing GO term enrichment analysis. I am personally using [clusterProfile](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html). Here the two publications:

* [Yu et al 2012](https://www.liebertpub.com/doi/10.1089/omi.2011.0118)
* [Wu et al 2021](https://www.sciencedirect.com/science/article/pii/S2666675821000667?via%3Dihub)

These two resources give nice examples of how to use clusterProfiler: [https://yulab-smu.top/biomedical-knowledge-mining-book/index.html](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html) and [https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/]([https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/).

I detailed below how to run a GO term analysis using B73 NAM5 gene ID (Zm00000ebXXXXXXX)

## Install clusterProfiler

Get all these libraries clusterProfiler needs:

```{r}
# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
bioconductor_packages <- c("DOSE","clusterProfiler","enrichplot","AnnotationDbi","GO.db")
BiocManager::install(bioconductor_packages)

# Install CRAN packages
cran_packages <- c("tidyverse","readxl","xlsx")
for (count in 1:length(cran_packages)) install.packages(cran_packages[count])

# Cross fingers and hope everything was installed properly
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(xlsx)
library(readxl)
library(AnnotationDbi)
library(GO.db)
library(tidyverse)

```

## Create TERM2NAME and TERM2GENE R objects

clusterProfiler needs 3 things to work:
* TERM2NAME data frame, which contains two columns: GO term ID (GO:0009880) and GO description "embryonic pattern specification"
* A TERM2GENE data frame, which contains two columns: GO term ID (GO:0000001) and B73 NAM5 gene ID "Zm00001eb000800"
* vector containing the B73 NAM5 gene IDs to be processed into the clusterProfiler function "enrich"

I have already added the TERM2NAME and TERM2GENE R objects in this GitHub repository, so you only need to generate TERM2GENE file if you work with a different gaf file than the one generated in this GOMAP pipeline for B73 NAM5 annotation.

### Prepare TERM2NAME object

Get a summary of all 44,509 GO ID and their terms (http://geneontology.org/docs/GO-term-elements). These are available in the GOTERM object of the AnnotationDbi package.

> Every term has a human-readable term name — e.g. mitochondrion, glucose transmembrane transport, or amino acid binding — and a GO ID, a unique seven digit identifier prefixed by GO:, e.g. GO:0005739, GO:1904659, or GO:0016597

This step is needed to be performed only once and can then be used across any GO analysis and any species.

```{r}
library(AnnotationDbi)

# Put GOTermsAnnDbBimap object into a dataframe
df_GOTERM <- as.data.frame(GOTERM)

# Remove first column (otherwise, two redundant column)
df_GOTERM <- df_GOTERM[,2:7]

# Create a list of dataframe with terms classified by ontology types
TERM2NAME_all <- df_GOTERM %>% dplyr::select(go_id,Term) %>% unique()
TERM2NAME_BP <- df_GOTERM %>% dplyr::filter(Ontology == "BP") %>% dplyr::select(go_id,Term) %>% unique()
TERM2NAME_CC <- df_GOTERM %>% dplyr::filter(Ontology == "CC") %>% dplyr::select(go_id,Term) %>% unique()
TERM2NAME_MF <- df_GOTERM %>% dplyr::filter(Ontology == "MF") %>% dplyr::select(go_id,Term) %>% unique()

# Create a list of dataframe
TERM2NAME <- list(ALL=TERM2NAME_all, BP=TERM2NAME_BP, CC=TERM2NAME_CC, MF=TERM2NAME_MF)

# Save this R object
saveRDS(TERM2NAME, "TERM2NAME.rds")

```

### Prepare TERM2GENE object

We now want to import the modified output gaf file generated by GOMAP to extract the geneID and GO terms.

```{r}
# Import data in R
TERM2GENE <- read.delim("maize_B73_NAM5.aggregate.GO_gene.gaf", header=F)
colnames(TERM2GENE) <- c("GO","gene")

saveRDS(TERM2GENE, "TERM2GENE.rds")
```

That's it!

Now we have everything set to perform our first GO term analysis.

# Run GO term enrichment analysis

I created a [go_functions.R] file that can be sourced in your current R or Rmarkdown file to call the functions needed to run the analysis.

The [go_functions.R] file contains two R functions:

The function `ego_analysis` performs three GO term analyses, one for each ontology: Biological process (BP), cellular component (CC), and molecular function (MF). One can explicit the wanted ontology when using the function `enrichGO` using the argument "ont" (e.g. `ont="MF"`) but this is not possible when using home-made GO annotation (as we do here for B73 NAM5). In the case of home-made annotation, clusterProfiler provides the function `enrich` but there is no possibility to define what ontology to investigate. I, therefore, split the TERM2NAME object into three separate objects, containing BP, CC, and MF GO terms. The home-made function `ego_analysis` will return a list of ego result objects (class enrichResult) that can be processed further for plotting with for instance the `dotplot` function of the `enrichplot` package (we will see that below). One can also easily turn these "ego" objects into data frames to perform manipulation with `tidyverse` functions for instance.

The function `go_search` is a practical tool to find correspondences between gene ID, GO ID, and a GO term name. Little examples are worth lengthy explanations:

```
# Get me all GO ID matching the maize gene Zm00001eb000690
head(go_search(method="gene2GO", "Zm00001eb000690"))
       go_id                                                          Term
1 GO:0000976 transcription regulatory region sequence-specific DNA binding
2 GO:0003700                     DNA-binding transcription factor activity
3 GO:0005634                                                       nucleus
4 GO:0006351                                  transcription, DNA-templated
5 GO:0006355                    regulation of transcription, DNA-templated
6 GO:0008289                                                 lipid binding

# Get me all maize genes being annotated with the GO ID GO:0003700
head(go_search(method="GO2gene","GO:0003700"))
          GO            gene
1 GO:0003700 Zm00001eb000210
2 GO:0003700 Zm00001eb000470
3 GO:0003700 Zm00001eb000690
4 GO:0003700 Zm00001eb001210
5 GO:0003700 Zm00001eb001460
6 GO:0003700 Zm00001eb001550

# Get the GO ID for the term description "nucleus"
go_search(method="term2GO", "nucleus")
       go_id    Term
1 GO:0005634 nucleus
```

Now, let's run our first analysis. You need this GitHub repository and all libraries specified above installed.

```{r}
# Clone this repository
git clone https://github.com/johanzi/GOMAP_maize_B73_NAM5.git

# Be sure to source the go_functions.R file properly
source("/path/to/GOMAP_maize_B73_NAM5go_functions.R")

# Create a list of 100 random genes
TERM2GENE <-  readRDS("TERM2GENE.rds")
random_100_geneID <- sample(TERM2GENE$gene, size=100)

# Run the enrichment analysis
# Let's select all genes annotated with
# GO ID GO:0019684 (photosynthesis, light reaction)
# That makes a vector of 85 genes
geneID_GO0019684 <- GO_analysis_data %>% filter(GO=="GO:0019684") %>% dplyr::select(gene) %>% pull()

# Run the ego_analysis function providing the vector of geneID argument
ego <- ego_analysis(geneID_GO0019684)

# Check if any of the three analyses yielded GO terms
# with significant p-adjusted value.
lapply(ego_analysis, function(x) sum(x@result$p.adjust < 0.05))

$ego_BP
[1] 112

$ego_CC
[1] 19

$ego_MF
[1] 14

# We get significant hits in all three ontologies
```
If any significant enrichment is found, one can visualize in GO term plot using `dotplot` function from the enrichplot package.

```{r}
dotplot(ego_analysis$ego_BP,showCategory=10, title="BP")
dotplot(ego_analysis$ego_CC,showCategory=10, title="CC")
dotplot(ego_analysis$ego_MF,showCategory=10, title="MF")
```

![](images/GO_analysis_plots.JPG)

Obviously, the term "GO:0019684 (photosynthesis, light reaction)" is the strongest signal since we selected all 85 genes annotated with this term.

To dig further into the GO enrichment analysis of BP ontology, one can turn it into a data frame:

```{r}
df_ego_analysis_BP <- as.data.frame(ego_analysis$ego_BP@result)

ego$ego_BP@result
```

Here it is, the end of the pipeline. I hope this will save time for others that struggle finding a GO term annotation and a straight forward way to visualize GO term enrichment analysis results.

## Authors

* **Johan Zicola** - [johanzi](https://github.com/johanzi)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

