FROM cyversevice/rstudio-verse:4.0.0-ubuntu18.04
# Pulls rocker tidyverse container with ubuntu 18.04 LTS for base build

# Build-time metadata as defined at http://label-schema.org
ARG BUILD_DATE
ARG VCS_REF
ARG VERSION
LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="RStudio Verse" \
      org.label-schema.description="Built from Rocker-Project RStudio Verse, additional depends for CyVerse K8s workbench" \
      org.label-schema.url="https://cyverse.org" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="e.g. https://github.com/cyverse-vice/rstudio-verse" \
      org.label-schema.vendor="CyVerse" \
      org.label-schema.version=$VERSION \
      org.label-schema.schema-version="1.0.0"

LABEL maintainer="Ryan Bartelme rbartelme@arizona.edu"
# This image builds off the rocker tidyverse image of rstudio
# and adds: cutadapt for primer trimming fastq files,
# as well as NEON api access to download fastq files,
# and dada2 for generating ASV's from fastqs

# ubuntu setup
RUN apt-get update && apt-get install -y dialog apt-utils tzdata

# install additonal R libraries
USER root
RUN R -e 'install.packages(c("BiocManager","neonUtilities", "vegan", "zetadiv", \
          "deblur", "decontam", "R.utils", "testthat"), repos = "https://cloud.r-project.org/", \
          dependencies = TRUE)'

# set BiocManager Version
RUN R -e 'BiocManager::install(version = "3.11", ask = FALSE)'

# install rhdf5 for neonUtilities
RUN R -e 'BiocManager::install("rhdf5")'
# install dada2 with biocmanager version 3.11
RUN R -e 'BiocManager::install("dada2", version = "3.11")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("DECIPHER")'
RUN R -e 'BiocManager::install("phyloseq")'

# install geoNEON
RUN R -e 'devtools::install_github("NEONScience/NEON-geolocation/geoNEON", dependencies=TRUE)'

#install neonMicrobe
RUN R -e 'devtools::install_github("claraqin/neonMicrobe", dependencies=TRUE)'

# setup rstudio user for conda/python/cutadapt
USER rstudio
WORKDIR /home/rstudio

# install miniconda
ENV PATH="/home/rstudio/miniconda3/bin:${PATH}"
ARG PATH="/home/rstudio/miniconda3/bin:${PATH}"
RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /home/rstudio/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# update miniconda3
RUN conda init bash \
  && . ~/.bashrc \
  && conda update -n base -c defaults conda \
  && conda config --add channels defaults \
  && conda config --add channels bioconda \
  && conda config --add channels conda-forge \
  && conda install -y cutadapt \
  && echo 'alias cutadapt = "conda run --prefix '/home/rstudio/miniconda3/envs/cutadaptenv' cutadapt"'

COPY profile .profile

USER root
ENTRYPOINT ["/usr/local/bin/run.sh"]
