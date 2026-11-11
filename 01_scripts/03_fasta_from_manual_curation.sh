# This script is here to obtain a mitochondrial sequence from a getOrganelle genome graph after manual curation in Bandage

#user variable:
GFA=$1
OUT=$2
NCPU=$3

if [[ -z "$NCPU" ]]
then
    NCPU=1
fi

get_organelle_from_reads.py -F animal_mt -g GFA -o $OUT 
