#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(Biostrings)
library(tidyverse)

# 16S and ITS were sequenced using illumina MiSeq 2x300
#===============================================================================
#  Set-Up
#===============================================================================

# detect number of logical cores, requires parallel R library
ncores <- detectCores(all.tests = TRUE, logical = TRUE)

# set logical cores to use in multi-thread based on available
if(ncores >= 48){
  core_use <- 48
}else{
  core_use <- ncores
}

#===============================================================================
# Begin ITS Amplicon Analysis of NEON Data
#===============================================================================

#set relative path for its fastq files
its_path <- "~/fastq/its"

#list full paths for forward and reverse reads, sorting alphanumerically 
fnFs_its <- sort(list.files(its_path, pattern = "_R1.fastq", full.names = TRUE))
fnRs_its <- sort(list.files(its_path, pattern = "_R2.fastq", full.names = TRUE))

#input forward and reverse ITS primers from NEON database
its_fwd <- "CTTGGTCATTTAGAGGAAGTAA" #forward
its_rev <- "GCTGCGTTCTTCATCGATGC" #reverse

#define function to look for primers in all orientations
all_orients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # Biostrings uses DNAString obj vs. character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

#find all orientations for each ITS primer
its_fwd_orients <- all_orients(its_fwd)
its_rev_orients <- all_orients(its_rev)

#check for filtN directory, if not present create it
filt_n <- file.path(its_path, "filtN")
if(!dir.exists(filt_n)) dir.create(filt_n)

#prefilter all ITS sequences with N's before primer trimming to speed up calcs
fnFs_its_filtn <- file.path(its_path, "filtN", basename(fnFs_its))
fnRs_its_filtn <- file.path(its_path, "filtN", basename(fnRs_its))
filterAndTrim(fnFs_its, fnFs_its_filtn, fnRs_its, fnRs_its_filtn,
                maxN = 0, multithread = TRUE)
# error due to untar not having finished in the dir
system('cd ~/fastq/its && for i in *.tar; do tar -xvf $i; done')
# files pre-processed but error thrown, therefore debug with setdiff list files
# Error message was...```Error in add(bin) : record does not start with '@'```

