#!/usr/bin/env bash

# Modified from Mike Lee's cutadapt example:
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex

# create a list of the files to trim pass to samples bash variable
ls *_R1.fastq | sed 's/_R1.fastq//' > samples

#prompt user for inputs
read -p 'Input Forward Primer: ' fwd
read -p 'Input Reverse Primer:' rev
read -p 'Specify cutadapt diagnostic output file name:' cud_out

# loop to trim primers
for sample in $(cat samples)
do

 echo "On sample: $sample"
    cutadapt -a $fwd \
    -A $rev \
    -m 215 -M 285 --discard-untrimmed \
    -o ${sample}_R1_trimmed.fastq -p ${sample}_R1_trimmed.fastq \
    ${sample}_R1.fastq ${sample}_R2.fq \
    >> ${cud_out}.txt 2>&1
done