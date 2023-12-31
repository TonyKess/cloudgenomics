# Mamba and Bioconda
This section details how to build the conda environments, download key software, and access supplementary information that is required for carrying out genomic analyses in a cloud environment.
This section mirrors much of the information provided in doing this work in an [HPC environment](https://github.com/TonyKess/genotyping_hpc/edit/main/env_building.md), the advantadge of the cloud environment being that we can work
directly with our software and data without working in a shared environment with compute-intensive analyses run using a shceduling system. 

[Mamba](https://github.com/mamba-org/mamba) is a fast and parallel implementation of the package management tool [conda](https://docs.conda.io/projects/conda/en/stable/#)
More information on installing conda can be found [here](https://github.com/GRDI-GenARCC/tutorials-and-workshops/blob/main/Conda/conda_installation_guide.md)

Mamba is a package manager that allows you to install and manage software packages. Many of the most frequently used bioinformatics tools are packaged and easily isntallable via [Bioconda](https://bioconda.github.io/).To use mamba to install software via Bioconda, you need to first install mamba.

You can install mamba using the following command:

```bash
conda install mamba -c conda-forge
```

Alternatively, you can download from conda-forge directly:

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
```

Now we activate and install:
```bash
chmod +x Mambaforge-Linux-x86_64.sh
./Mambaforge-Linux-x86_64.sh
```

During installation we can select install location - on a cloud VM it's probably best to use one of the installed datadisks to ensure there's plenty of space available
```
/datadisk0/software/mambaforge/
```

# Building a genotyping environment

Now we can start installing software needed! 
Ths first thing we'll want to do with our data is run quality and adapter trimming. We can use fastp to do this

```
mamba create -n fastp fastp -c bioconda
```

We will also make an environment for alignment of reads to the reference genome, as well as working with genotype likelihoods in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD). These are imported together to ensure version comaptibility between ANGSD and the corresponding htslib and samtools installs. 

```
mamba create -n align samtools htslib bwa bwa-mem2 angsd -c bioconda
```
To get the current version of GATK, we do:

```
mamba create -n gatk4 gatk4 -c bioconda
```

Next we need different versions of software for genotyping and pre-genotyping processes, that are not present in useful versions in bioconda. I like to store these alongside mamba environments in:
```
/software
```

ANGSD works with inputs from older versions of GATK, whereas we can directly carry out genotyping in the current version of GATK. To get the old verison (3.7), we do:

```
wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2 &&
bzip2 -d GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2 &&
tar -xvf GenomeAnalysisTK-3.7-0-gcfedb67.tar
```
This version of GATK also needs its own, older version of java that we can install using conda. 
```
mamba create -n gatk37 openjdk=8
```
We can now move on to a genotyping workflow that we've used in the past in HPC environments. 
