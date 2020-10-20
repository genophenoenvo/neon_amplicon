# NEON Microbiome Amplicon Analyses

---

## Software Included in Container Environment

The Docker container is built from the R 4.0 Ubuntu 18.04 Tidyverse Rstudio Server [base image](https://hub.docker.com/layers/rocker/tidyverse/4.0.0-ubuntu18.04/images/sha256-f8947c19fda376764dcc35fcb201e5567826474ba4bfcde54bf5424e2225fa2f?context=explore) created by the [Rocker Group](https://www.rocker-project.org/). 

This analysis generally follows Mike Lee's suggested [dada2](https://benjjneb.github.io/dada2/) [workflow](https://astrobiomike.github.io/amplicon/dada2_workflow_ex). Including the use of [cutadapt](https://cutadapt.readthedocs.io/) for primer trimming. As per Mike's tutorial, `cutadapt` is run from the terminal, with the notable exception that a shell script is provided to iterate over the fastq files pulled from the NEON database. 

---

## Docker Container Rstudio Environment Setup

### Pulling the container from dockerhub

The Docker container for these analyses can be run on CyVerse VICE as an application, or you may pull the container onto any compatible machine with:

`docker pull cyversevice/rstudio-neon-dada2:1.1`

### Running the container image



### Logging in to Rstudio Server



### Initial Setup

There is some minor house keeping to do in the Rstudio terminal tab before pulling this repository into the container environment.

#### 1) Set up iRODS


#### 2) Setup .Renviron file with NEON API Key

### Pulling this repository into the container

---

## Running the Analyses in this Repository

