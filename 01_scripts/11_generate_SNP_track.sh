#!/bin/bash


# Generate SNP count in windows size from pileup files for each species in pileup

wins=500
pileup="12_pileup_files"
out="13_tracks"

for i in $pileup/*mpileup
do
python3 01_scripts/00_individual_utils/create_SNP_density_track.py $i $out 500
done
