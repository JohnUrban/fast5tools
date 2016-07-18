#!/bin/bash

if [ $# -eq 0 ]; then echo "
Usage:
fast5DerivedFastaReadStats.sh file.fa

Output - tab-separated:
Molecule_name, read_type, read_length, mean_quality_score
"; exit; fi

grep ">" $1 | awk 'OFS="\t" {sub(/>/,""); gsub(/\|/,"\t"); gsub(/:/,"\t"); print $11"_"$13"_"$7"_"$9, $1, $3, $5}'
