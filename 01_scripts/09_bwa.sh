#!/bin/bash

# Align raw reads to mitochondrial sequences
NCPU=$1

reads="02_raw_files"
bwa_out="11_bam_alignements"
fasta_in="06_reoriented"
misc_paper="99_paper_specific_files/Stroggylo_extra"
SNPs="12_pileup_files"

cat $fasta_in/Diploconelus_abdominalis.fasta $fasta_in/Tomosvaryella_diploconelus_abdominalis.fasta > $bwa_out/Diploconelus_abdominalis.fasta
cat $fasta_in/Psammotettix_alienus.fasta $fasta_in/Sorhoanus_sp.fasta > $bwa_out/Psammotettix_alienus.fasta
cat $fasta_in/Stroggylocephalus_mixtus.fasta $misc_paper > $bwa_out/Stroggylocephalus_mixtus.fasta
cp $fasta_i Hebecephalus_algidus.fasta $bwa_out/
cp $fasta_i Limotettix_myrialis.fasta $bwa_out/
cp $fasta_i Limotettix_plutonius.fasta $bwa_out/

for i in $bwa_out/*fasta
  do
  sp=$(basename ${i/.fasta})
  bwa_index $i
  bwa mem -t NCPU $i $reads/sp_R1.fq.gz $reads/sp_R2.fq.gz | \
  samtools sort - \
  samtools view -q 10 -b -o $bwa_out/${sp}.bam -

  bcftools mpileup -Ou -f $i -o $SNPs/${sp}.mpileup $bwa_out/${sp}.bam
  rm $i
  done
