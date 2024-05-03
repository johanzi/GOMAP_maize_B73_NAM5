# GOMAP_maize_B73_NAM5
GO term annotation of the maize B73 NAM5 assembly

# Objectives

I needed to perform GO enrichment analysis on genes targeted by sRNAs in maize before the time a B73 NAM5 GO term annotation was available on [MaizeGDB](https://download.maizegdb.org/GeneFunction_and_Expression/Pannzer_GO_Terms/) in 2022. This annotation was performed using [PANNZER](http://ekhidna2.biocenter.helsinki.fi/sanspanz/), see reference paper [Törönen & Holm, 2022](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4193). I therefore decided to use the Gene Ontology Meta Annotator for Plants ([GOMAP]([https://bioinformapping.com/gomap/master/RUNNING.htm](https://bioinformapping.com/gomap/master/OVERVIEW.html)) pipeline.

Both PANNZER and GOMAP annotated 39,756 genes. However, GOMAP could generate 493,310 annotations vs 167,519 for PANNZER. Although GOMAP manual recommends running the pipeline on a high-performance cluster (HPC) due to the sheer amount of calculation it needs, the setting of my HPC provider did not allow me to run such a long job. It took my local workstation (10 cores of 2,2 GHz each) about one month to complete the annotation. To save time and energy for others, I provide here the GOMAP annotation for B73 NAM5 with its step-by-step.

# Requisites

* UNIX-based computer with a decent amount of RAM and cores (mine had 10 x v4 processors of 2.2 GHz each, 63 Gb of RAM, 4 Gb GPU, OS Linux Mint 19.1 Cinnamon)
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
* Your email, needed to run one of the tool (Argot GO)
* cpus, depends on your workstation, the more, the faster your job will be. Plan anyway a good month to let your machine run if you don't have access to a HPC



