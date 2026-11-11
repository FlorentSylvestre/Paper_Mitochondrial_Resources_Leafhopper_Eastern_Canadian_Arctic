!#/bin/bash

MITOSGFF="06_reoriented"
TRNAGFF="07_tRNAscan_SE"

for i in $MITOSGFF/*gff
do
  spes=$(basename ${i/.gff})
  01_scripts/00_individual_utils/merge_mitos_tRNA.py $i $TRAGFF/${spes}.gff
done
