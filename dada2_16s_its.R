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
filt_out <- filterAndTrim(fnFs_its, fnFs_its_filtn, fnRs_its, fnRs_its_filtn,
                maxN = 0, multithread = core_use, verbose = TRUE)

# list files in its directory and compare those to the ones in the filtered dir
raw_its_fastq <- list.files(path = its_path, pattern = "fastq")
nfiltered_its_fastq <- list.files(path = file.path(its_path, "filtN"), 
                                  pattern = "fastq")
bad_header_fastq <- setdiff(raw_its_fastq, nfiltered_its_fastq)

# 6 files will not process due to line errors in the sequencing run, these are
# most likely corrupted files and should not be processed further
# offending files are stored in the variable bad_header_fastq & the rest of the
# files will be processed further

#check distribution of the orientation of primers in the first sample
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(its_fwd_orients, primerHits, fn = fnFs_its_filtn[[20]]), 
      FWD.ReverseReads = sapply(its_fwd_orients, primerHits, fn = fnRs_its_filtn[[20]]), 
      REV.ForwardReads = sapply(its_rev_orients, primerHits, fn = fnFs_its_filtn[[20]]), 
      REV.ReverseReads = sapply(its_rev_orients, primerHits, fn = fnRs_its_filtn[[20]]))

#cutadapt sanity check
system2(command = "cutadapt", args = "--version")

# take n filtered data and trim ITS primers

#create directory for cutadapt output
path_cut <- file.path(filt_n, "cutadapt")
if(!dir.exists(path_cut)) dir.create(path_cut)

#prepare file outputs
fnFs_its_cut <- file.path(path_cut, basename(fnFs_its_filtn))
fnRs_its_cut <- file.path(path_cut, basename(fnRs_its_filtn))

#store reverse compliments of primers as strings
its_fwd_rc <- dada2::rc(its_fwd)
its_rev_rc <- dada2::rc(its_rev)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
r1_flags <- paste("-g", its_fwd, "-a", its_rev_rc) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
r2_flags <- paste("-G", its_rev, "-A", its_fwd_rc) 

for(i in seq_along(fnFs_its_filtn)){
  system2(command = "cutadapt", args = c(r1_flags, r2_flags, "-n", 2, 
                    # -n2 required to remove both primers from reads
                    "-j", core_use, #multithread with cores detected in first lines
                    "--discard-untrimmed", #discard untrimmed reads
                    "-o", fnFs_its_cut[i], "-p", fnRs_its_cut[i], # output files
                    fnFs_its_filtn[i], fnRs_its_filtn[i])) # input files
}

# get paths for cutadapt trimmed files
cutFs_its <- sort(list.files(path_cut, pattern = "_R1.fastq", full.names = TRUE))
cutRs_its <- sort(list.files(path_cut, pattern = "_R2.fastq", full.names = TRUE))

#sanity check for primer cutting
rbind(FWDp_ForwardReads = sapply(its_fwd_orients, primerHits, fn = cutFs_its[[20]]), 
      FWDp_ReverseReads = sapply(its_fwd_orients, primerHits, fn = cutRs_its[[20]]), 
      REVp_ForwardReads = sapply(its_rev_orients, primerHits, fn = cutFs_its[[20]]), 
      REVp_ReverseReads = sapply(its_rev_orients, primerHits, fn = cutRs_its[[20]]))



#extract sample name from fwd read file names
get_sample_name <- function(fname) strsplit(basename(fname), "_R1.fastq")[[1]][1]
its_original_sample_names <- unname(sapply(cutFs_its, get_sample_name))

#sanity check
head(its_sample_names)
head(cutFs_its)

#check read quality before filtering and trimming
plotQualityProfile(cutFs_its[1:2])
plotQualityProfile(cutRs_its[1:2])

#filter and trim output directory creation
its_filt_path <- file.path(path_cut, "filtered")
if(!dir.exists(its_filt_path)) dir.create(its_filt_path)

#filtered file output assignments
filtFs_its <- file.path(its_filt_path, basename(cutFs_its))
filtRs_its <- file.path(its_filt_path, basename(cutRs_its))

#filter and trim reads
its_filt_out <- filterAndTrim(cutFs_its, filtFs_its, cutRs_its, filtRs_its, 
                     maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50,
                     rm.phix = TRUE, compress = TRUE, multithread = core_use)

#sanity check output
head(its_filt_out)

#since it's 2x300bp chemistry, and multiplexed, 
  # it's not surprising to lose 50% of the reads

#learn forward error rates
errF_its <- learnErrors(filtFs_its, multithread = core_use)

#learn reverse error rates
errR_its <- learnErrors(filtRs_its, multithread = core_use)

# need to reset filtered file names & paths,
# since 3 samples were dropped due to quality
filtered_fwd_its <- list.files(path = its_filt_path,
                               pattern = "R1.fastq", full.names = TRUE)
filtered_rev_its <- list.files(path = its_filt_path,
                               pattern = "R2.fastq", full.names = TRUE)

#dereplicate reads
derepFs_its <- derepFastq(filtered_fwd_its, verbose = TRUE)
derepRs_its <- derepFastq(filtered_rev_its, verbose = TRUE)
# Name the derep-class objects by the sample names
# replace old sample character vector
its_sample_names <- unname(sapply(filtered_fwd_its, get_sample_name))

#assign sample names to dereplicated reads
names(derepFs_its) <- its_sample_names
names(derepRs_its) <- its_sample_names

#sample inference
dadaFs_its <- dada(derepFs_its, err = errF_its, multithread = core_use)
dadaRs_its <- dada(derepRs_its, err = errR_its, multithread = core_use)

#merged paired reads
merged_its <- mergePairs(dadaFs_its, derepFs_its, dadaRs_its,
                            derepRs_its, verbose=TRUE)

#make sequencing table
seqtab_its <- makeSequenceTable(merged_its)
dim(seqtab_its)

#remove chimeras
nochimera_seqtab_its <- removeBimeraDenovo(seqtab_its, method = "consensus",
                                           multithread = core_use, verbose = TRUE)


dim(seqtab_its)
distribution_its <- as.data.frame(table(nchar(getSequences(nochimera_seqtab_its))))
colnames(distribution_its) <- c("read_length", "counts")
barplot(height = distribution_its$counts, names.arg = distribution_its$read_length,
        xlab = "Merged Read Length", ylab = "Frequency", )
p<-ggplot(data=distribution_its, aes(x=read_length, y=counts)) +
  geom_bar(stat="identity")
p+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#track reads through pipeline 
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs_its, getN), sapply(dadaRs_its, getN),
               sapply(merged_its, getN), rowSums(nochimera_seqtab_its))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- its_sample_names
head(track)
tracked <- as.data.frame(cbind(track, track[,4]/track[,1]))
colnames(tracked)[5] <- "percent_filtered_reads"
tracked

#===============================================================================
# Begin 16S Amplicon Analysis of NEON Data
#===============================================================================

#set relative path for its fastq files
path_16s <- "~/fastq/16s"

#list full paths for forward and reverse reads, sorting alphanumerically 
fnFs_16s <- sort(list.files(path_16s, pattern = "_R1.fastq", full.names = TRUE))
fnRs_16s <- sort(list.files(path_16s, pattern = "_R2.fastq", full.names = TRUE))

#input forward and reverse 16S primers from NEON ***need to update values to V3-V4
fwd_16s <- "CTTGGTCATTTAGAGGAAGTAA" #forward
rev_16s <- "GCTGCGTTCTTCATCGATGC" #reverse

#find all orientations for each 16S primer
fwd_16s_orients <- all_orients(fwd_16s)
rev_16s_orients <- all_orients(rev_16s)

#create output directory for cutadapt
cut_16s <- file.path(path_16s, "cutadapt")
if(!dir.exists(cut_16s)) dir.create(cut_16s)




