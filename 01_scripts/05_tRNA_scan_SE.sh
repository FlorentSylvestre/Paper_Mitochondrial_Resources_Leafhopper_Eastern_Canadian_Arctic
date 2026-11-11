#!/bin/bash

# run tRNA_scan_SE on every mitochondria to complete tRNA annotation
# THIS CODE ASSUME THAT TRNASCAN IS INSTALLED IN A CONDA ENV. EDIT THE GCODE PATH TO MATCH YOUR SYSTEM!!!

GCODE=$CONDA_PREFIX/share/tRNAscan-SE/gcode/gcode.invmito

for i in 06_reoriented/*fasta
  do
  tRNA_scan_SE -O -o 07_tRNAscan_SE/$i{./fasta}.out -j 07_tRNAscan_SE/${i/.fasta}.gff -g $GCODE $i
  done
