#!/bin/bash

#Usage: 02_getOrganelle.sh NCPU
#This software is best install in a conda environment:
# conda create -n getOrganelle
# conda activate getOrganelle
# conda install conda install bioconda::getorganelle-1.7.7.1

FASTP_PATH="02_fastp"
OUT_DIR="03_getOrganelle"

#user variable:
NCPU=$1
if [[ -z "$NCPU" ]]
then
    NCPU=1
fi

ls $"FASTP_PATH"*_1.fastq.gz | sed "s/./\\///" sed "s/*_1.fastq.gz//"| \
  parallel -j $NCPU get_organelle_from_reads.py -1 {}_1.fastq.gz \
  -2 {}_2.fastq.gz -F animal_mt -t 1 -o $"OUT_DIR"/{}/ \
  -k 21,45,65,85,105

