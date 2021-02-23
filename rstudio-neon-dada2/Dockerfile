FROM rocker/verse:4.0.0-ubuntu18.04

# Pulls rocker tidyverse container with ubuntu 18.04 LTS for base build

MAINTAINER "Ryan Bartelme rbartelme@arizona.edu"
# This image builds off the rocker tidyverse image of rstudio
# and adds: cutadapt for primer trimming fastq files,,
# as well as NEON api access to download fastq files,
# and dada2 for generating ASV's from fastqs

## Install CyVerse VICE Depends

RUN apt-get update && apt-get install -y lsb-base lsb-release wget \
    apt-transport-https python2.7 python-requests curl supervisor \
    nginx gnupg2 dialog apt-utils tzdata

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
RUN R -e 'install.packages(c("BiocManager","neonUtilities", "vegan", "zetadiv", \
          "deblur", "decontam", "R.utils"), repos = "https://cloud.r-project.org/", \
          dependencies = TRUE)'

#set BiocManager Version
RUN R -e 'BiocManager::install(version = "3.11", ask = FALSE)'

#install rhdf5 for neonUtilities
RUN R -e 'BiocManager::install("rhdf5")'
# install dada2 with biocmanager version 3.11
RUN R -e 'BiocManager::install("dada2", version = "3.11")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("DECIPHER")'
RUN R -e 'BiocManager::install("phyloseq")'

#install geoNEON
RUN R -e 'devtools::install_github("NEONScience/NEON-geolocation/geoNEON", dependencies=TRUE)'

# install python 3 setup utilities pip3 and tzdata
RUN apt-get install -y dialog apt-utils tzdata

#setup rstudio user for conda/python/cutadapt
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
  && conda install cutadapt \
  && echo 'alias cutadapt = "conda run --prefix '/home/rstudio/miniconda3/envs/cutadaptenv' cutadapt"'

#irods path creation
RUN mkdir /home/rstudio/.irods/

# irods json config file
RUN echo '{"irods_host": "data.cyverse.org", "irods_port": 1247, "irods_user_name": "", "irods_zone_name": "iplant"}' >> /home/rstudio/.irods/irods_environment.json

#------------------------#
# Configs for deployment #
#------------------------#

ENV PASSWORD "rstudio1"
RUN bash /etc/cont-init.d/userconf

COPY run.sh /usr/local/bin/run.sh
RUN chmod +x /usr/local/bin/run.sh

COPY nginx.conf.tmpl /nginx.conf.tmpl
COPY rserver.conf /etc/rstudio/rserver.conf
COPY supervisor-nginx.conf /etc/supervisor/conf.d/nginx.conf
COPY supervisor-rstudio.conf /etc/supervisor/conf.d/rstudio.conf

ENV REDIRECT_URL "http://localhost/"

# add privileges to rstudio user
ARG LOCAL_USER=rstudio
ARG PRIV_CMDS='/bin/ch*,/bin/cat,/bin/gunzip,/bin/tar,/bin/mkdir,/bin/ps,/bin/mv,/bin/cp,/usr/bin/apt*,/usr/bin/pip*,/bin/yum'

RUN if [ -x /usr/bin/apt ]; then \
      apt-get update && apt-get -y install sudo && rm -rf /var/lib/apt/lists/*; \
    elif [ -x /bin/yum ]; then \
      yum -y update && yum -y install sudo && yum clean all; \
    fi

#eliminate password at login
RUN echo "$LOCAL_USER	ALL=NOPASSWD: $PRIV_CMDS" >> /etc/sudoers

#change appearance of shell prompt
USER rstudio
RUN echo 'export PS1="[\u@cyverse] \w $ "' >> /home/rstudio/.bash_profile
USER root

#entrypoint
ENTRYPOINT ["/usr/local/bin/run.sh"]
