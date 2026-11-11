#!/bin/bash

# Script that shift GFF coordinates, then sort and prepare GFF for final submission to NCBI and analysis

out="10_final_gff"
coord_shifts="99_paper_specific_files/start_stop_changes.txt"
gffin="08_merged_annotations"

for i in gffin/*gff
    do
    python3 01_scripts/00_individual_utils/Shift_start_stop_codon.py $i $coord_shifts $out/$(basename $i)
    python3 01_scripts/00_individual_utils/Sort_clean_gff.py $coord_shifts $out/$(basename $i)
    done
