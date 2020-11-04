#!/usr/bin/env bash

#decompress all gzip in all subdirectories
cd ~/fastq && find . -name '*.gz' -execdir gzip -d '{}' \;
