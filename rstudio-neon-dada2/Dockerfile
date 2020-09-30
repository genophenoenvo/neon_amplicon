FROM rocker/verse:4.0.0-ubuntu18.04

# Pulls rocker tidyverse container with ubuntu 18.04 LTS for base build

MAINTAINER "Ryan Bartelme rbartelme@arizona.edu"
# This image builds off the rocker tidyverse image of rstudio
# and adds: cutadapt for primer trimming fastq files,,
# as well as NEON api access to download fastq files,
# and dada2 for generating ASV's from fastqs

## Install CyVerse VICE Depends
## lsb is deprecated
RUN apt-get update && apt-get install -y lsb-base lsb-release wget apt-transport-https python2.7 python-requests curl supervisor nginx gnupg2

RUN apt-get install libfuse2

RUN curl "http://ftp.se.debian.org/debian/pool/main/o/openssl/libssl1.0.0_1.0.1t-1+deb8u8_amd64.deb" -O && \
    dpkg -i libssl1.0.0_1.0.1t-1+deb8u8_amd64.deb && \
    rm libssl1.0.0_1.0.1t-1+deb8u8_amd64.deb

# Install iCommands
RUN wget https://files.renci.org/pub/irods/releases/4.1.10/ubuntu14/irods-icommands-4.1.10-ubuntu14-x86_64.deb && dpkg -i irods-icommands-4.1.10-ubuntu14-x86_64.deb

ADD https://github.com/hairyhenderson/gomplate/releases/download/v2.5.0/gomplate_linux-amd64 /usr/bin/gomplate
RUN chmod a+x /usr/bin/gomplate

# provide read and write access to Rstudio users for default R library location
RUN chmod -R 777 /usr/local/lib/R/site-library

# install additonal R libraries
RUN R -e 'install.packages(c("BiocManager","neonUtilities", "vegan", "zetadiv", "deblur", "decontam"), dependencies = TRUE)'

# install dada2 with biocmanager version 3.11
RUN R -e 'BiocManager::install("dada2", version = "3.11")'

# install python 3 setup utilities pip3
RUN apt-get install -y python3-setuptools python3-pip

# install miniconda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# update miniconda3
RUN conda update -n base -c defaults conda
# add default channels
RUN conda config --add channels defaults
# add bioconda channels
RUN conda config --add channels bioconda
# add conda-forge channels
RUN conda config --add channels conda-forge
#install cutadapt via conda
RUN conda install cutadapt
#create cutadapt venv
RUN conda create -n cutadaptenv cutadapt

RUN chown -R rstudio /root/miniconda3
RUN chown -R rstudio /root/.conda

ENV PASSWORD "rstudio1"
RUN bash /etc/cont-init.d/userconf

COPY run.sh /usr/local/bin/run.sh
RUN chmod +x /usr/local/bin/run.sh

COPY nginx.conf.tmpl /nginx.conf.tmpl
COPY rserver.conf /etc/rstudio/rserver.conf
COPY supervisor-nginx.conf /etc/supervisor/conf.d/nginx.conf
COPY supervisor-rstudio.conf /etc/supervisor/conf.d/rstudio.conf

# iCommands json configuration for rstudio user profile
USER rstudio
RUN mkdir /home/rstudio/.irods/
RUN echo '{"irods_host": "data.cyverse.org", "irods_port": 1247, "irods_user_name": "", "irods_zone_name": "iplant"}' >> /home/rstudio/.irods/irods_environment.json
USER root


ENV REDIRECT_URL "http://localhost/"

ENTRYPOINT ["/usr/local/bin/run.sh"]
