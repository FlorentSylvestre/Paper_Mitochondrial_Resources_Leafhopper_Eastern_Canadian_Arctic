#!/bin/bash

new="15_new_gene_sequences"
acc="16_accession_gene_sequences"
concat="17_merged_gene_seq"
aligned="18_aligned_sequences"
# Merge outgroup and sample in one fasta for each gene
for i in $(ls $acc/)
do
        cat $new/$i $acc/$i > $concat/$i
        echo $i 
        mafft --thread 10 --maxiterate 1000 --localpair $concat/${i} >$aligned/$i
done
seqkit concat $aligned/* > $aligned/concatenated_alignement.fasta


