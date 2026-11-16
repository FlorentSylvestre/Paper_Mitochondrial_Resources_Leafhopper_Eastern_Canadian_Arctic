#!/bin/bash

# Usage: ./12_extract_prot_sequences.sh


acc_path="14_accessions"
acc_data="99_paper_specific_files/16_Accession_metadata.tar.gz"
acc_names="99_paper_specific_files/acc_taxtable"
new_seq_path="10_final_gff"
out_acc="16_accession_gene_sequences"
out_new="15_new_gene_sequences"

# Initialize files
tar -xvzf $acc_data -C $acc_path
cat $new_seq_path/*.fasta > tmp.fasta
cat $new_seq_path/*gff > tmp.gff
echo "Tomosvaryella_sp"> tmp.blacklist

python3 01_scripts/00_individual_utils/extract_prot_from_gff.py -f tmp.fasta -g tmp.gff -o $out_new -blacklist tmp.blacklist
rm  tmp.fasta tmp.gff tmp.blacklist

python3 01_scripts/00_individual_utils/extract_prot_from_gff.py -f $acc_path/Used_sequence.fasta -g $acc_path/Usedsequence.gff3 -o $out_acc -taxtable $acc_names



 
