#!/bin/python3

"""
Script that take a gff and a four column table with adjustement to start 
and stop positions for genes to calculate new position and output a new gff

example of annotation file (must be table separated): 
Diplocolenus_abdominalis    ATP6    15  0
Diplocolenus_abdominalis    Cox2    -15  -9
Diplocolenus_abdominalis    NAD3    -9  0
Diplocolenus_abdominalis    NAD4    30  62
Diplocolenus_abdominalis    NAD4L    -18  0

Usage: python3 adjust_start_stop_codon.py <path_to_gff> <path_to_Annotation_changes> <output_path>
"""

import sys
from collections import defaultdict

#Import args
gff_p, annotation_p, outfile_p = sys.argv[1:]

pos1 = 3
pos2 = 4
strand = 6


# parse annotation_p:
annotation_dict = defaultdict(lambda: defaultdict(list))

with open(annotation_p, "r") as handle:
    for line in handle:
        sp, gene, start, stop  = line.strip().split()
        annotation_dict[sp][gene] = [int(start), int(stop)]

# Parse and output gff
outfile = open(outfile_p, "w")
with open(gff_p, "r") as gffin:
    for line in gffin:
        if line.startswith("#"):
            outfile.write(line)
            continue

        l = line.strip().split()
        Name = [x for x in l[8].split(";") if x.startswith("Name")]
        Name = Name.pop().split("=")[1]
        
        if Name not in annotation_dict[l[0]]:
            outfile.write(line)
            continue

        if l[strand] == "-":
            start, stop = pos2, pos1 
        else:
            start, stop = pos1, pos2
        

        shift_s, shift_e = annotation_dict[l[0]][Name]
        l[start] = str(int(l[start]) + shift_s)
        l[stop] = str(int(l[stop]) + shift_e)
        outfile.write("\t".join(l) + "\n")

