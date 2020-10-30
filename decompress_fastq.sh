#!/usr/bin/env bash


# decompress all tar.gz in all subdirectories
cd ~/fastq && find . -name '*.tar.gz' -execdir tar -xzvf '{}' \;

#decompress all gzip in all subdirectories
cd ~/fastq && find . -name '*.gz' -execdir gzip -d '{}' \;
