#!/bin/bash

# Extract gff and fasta file flanking region for each gene to facilitate start and stop codon manual inspection

ls 08_merged_annotations/*gff >gff_list
ls 06_reoriented/*fasta >fasta_list
python3 01_scripts/00_individual_utils/Extract_flanking_regions.py fasta_list gff_list 500 09_manual_inspections
rm fasta_list gff_list
