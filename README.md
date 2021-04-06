# NEON Microbiome Amplicon Analyses

---

## Software Included in Container Environment

The Docker container is built from the R 4.0 Ubuntu 18.04 Tidyverse Rstudio Server [base image](https://hub.docker.com/layers/rocker/tidyverse/4.0.0-ubuntu18.04/images/sha256-f8947c19fda376764dcc35fcb201e5567826474ba4bfcde54bf5424e2225fa2f?context=explore) created by the [Rocker Group](https://www.rocker-project.org/). 

This analysis generally follows Mike Lee's suggested [dada2](https://benjjneb.github.io/dada2/) [workflow](https://astrobiomike.github.io/amplicon/dada2_workflow_ex). Including the use of [cutadapt](https://cutadapt.readthedocs.io/) for primer trimming. As per Mike's tutorial, `cutadapt` is run from the terminal, with the notable exception that a shell script is provided to iterate over the fastq files pulled from the NEON database. 

---

## Docker Container Rstudio Environment Setup

### Pulling the container from dockerhub

The Docker container for these analyses can be run on CyVerse VICE as an application, or you may pull the container image onto any compatible machine with:

`docker pull cyversevice/rstudio-neon-dada2:1.4`

### Running the container image



### Logging in to Rstudio Server



### Initial Setup

There is some minor house keeping to do in the Rstudio terminal tab before pulling this repository into the container environment.

#### 1) Set up iRODS

- **NOTE: You should have a CyVerse Account *before* running these commands**
- Navigate to the Terminal tab in the lower lefthand corner of RStudio.
- After selecting the Terminal tab type in `iinit` to configure iRODS
- When prompted for a username, enter your CyVerse user name
- And, similarly, when prompted for a password, enter your CyVerse password

This now gives access to the CyVerse data store associated with your account, and can allow for the use of [iCommands](https://learning.cyverse.org/projects/data_store_guide/en/latest/step2.html) to transfer data from this container to the CyVerse Data Store.

#### 2) Setup .Renviron file with NEON API Key

- **NOTE: You should have a NEON Account *before* generating the API Key Configuration**
- NEON API key generation tutorial can be found [here](https://www.neonscience.org/neon-api-tokens-tutorial)
- After getting your NEON API key it's time to create the `.Renviron` file
- Navigate to the terminal in Rstudio if you are not already there
- Run `vi .Renviron`
- This opens the vi text editor
- See the following [tutorial](http://heather.cs.ucdavis.edu/~matloff/UnixAndC/Editors/ViIntro.html) if you are unfamilar with vi operations, it's quite different than nano
- The following is an example of how you can write the API key into the environment:

```
NEON_TOKEN=WhateverYourGeneratedAPIkeyIS

```

- Do not quote your API key or use any spaces for the assignment to the NEON_TOKEN R environmental variable
- The empty line after the API key in the `.Renviron` file is important for this to work.
- Once this file is written, you will have full download speeds when running the scripts to query the NEON Data Portal.

### Pulling this repository into the container

Now this repository can be pulled into the container environment for a fully reproducible analysis.

- Navigate to File => New R Project
- Select "Version Control" from the Create Project window
- Then Select "Git" from the Create Project from Version Control window
- Fill in the fields as follows:

Repository url: `https://www.github.com/genophenoenvo/neon_amplicon/`
Project Directory: `neon_amplicon` (NOTE: This should automatically fill in after entering the url)
Create project as a subdirectory of: `~` (This should pull the repo into `/home/rstudio/`)

Now this repository should be inside the container and you can move to the Running the Analyses in this Repository Section of the README.

---

## Running the Analyses in this Repository

