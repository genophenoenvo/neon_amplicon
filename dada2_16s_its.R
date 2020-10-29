#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(Biostrings)

# 16S and ITS were sequenced using illumina MiSeq 2x300

#===============================================================================
# Begin ITS Amplicon Analysis of NEON Data
#===============================================================================

#set relative path for its fastq files
its_path <- "~/fastq/its"

#list full paths for forward and reverse reads, sorting alphanumerically 
fnFs_its <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs_its <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

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


