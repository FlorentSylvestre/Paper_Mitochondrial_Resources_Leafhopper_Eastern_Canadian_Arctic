#!/bin/bash
# Script that take a circular mitochondria fasta and associated gff
# And Pivot the coordinate to the start of selected gene.
# Script also Reverse complement sequence to put COX1 on forward strand if needed

for seq in 05_sequences_and_annotations*.fasta
  do
    name=${seq/.fasta}
    python3 01_scripts/00_individual_utils/reorient_MT.py seq ${name}.gff "trnI(gat)" 06_reoriented/${name}
  done

