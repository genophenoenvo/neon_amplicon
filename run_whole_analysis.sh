#!/usr/bin/env bash

# Prompt user for NEON API Key

# Write .Renviron file with NEON API Key

# run data pull, file hierarchy generation and clean-up
Rscript ~/neon_amplicon/neon_query/microbe_metadata.R

# Trim primers with cutadapt shell script

# Run dada2 analysis

# Analyze dada2 output
