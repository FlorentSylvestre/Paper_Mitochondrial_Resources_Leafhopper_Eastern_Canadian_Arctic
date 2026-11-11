#!/bin/bash

# Uses mosdepth to generate coverage track that can be parsed by proksee

output="13_tracks"
bam="12_pileup_files"
wins="500"
header="Start\tScore"

for i in $bam/*bam
  do
  prefix="$output/$(basename ${$i/.bam})"
  mosdepth_d4 -by $wins -m $prefix $i
  rm $prefix.mosdepth.global.dist.txt $prefix.mosdepth.region.dist.txt $prefix.mosdepth.summary.txt $prefix.per-base.bed.gz ${prefix}*.csi
  gunzip $prefix.region.bed.gz

  cut -f1,2,4 $prefix.region.bed >tmp
  for sp in $(cut -f1 tmp |sort | uniq)
  do
    echo -e $header >$output/$sp.covtrack
    grep -P "^$sp\t" tmp | cut -f 2,3 >> $sp.covtrack
  done
  rm $prefix.region.bed tmp
