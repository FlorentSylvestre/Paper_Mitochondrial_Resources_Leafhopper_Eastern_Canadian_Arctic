#!/bin/python

import sys
from collections import defaultdict

"""
Python script that take a fasta file, a gff annotation file from mitos2
and a gene name.

relocate the start position jsute before the gene name, and orient mitochondria 
so that COX 1 gene is on the forward strand.

Usage:
python3 reorient_MT.py <fasta_path>, <gff_path>, <pos>, <outprefix> 
fasta path: path of the fasta file to reorient
gff path: assosiated gff
pos: gene name OR physical position
outprefix: path to write new file

"""


class gff_entry:

    def __init__(self, line: None):
        self.line = line

        if line[0].startswith("#"):
            self.type = "header"
            self.description = {"Name":"header"}
            
        else:    
            self.type = line[2]
            self.start = int(line[3])
            self.end = int(line[4])
            self.strand = line[6]
            self.description = {x.split("=")[0]:x.split("=")[1] for x in line[8].split(";")}

    def get_line(self):
        if self == "header" or self == "region":
            return self.line 
        l = self.line
        l[3] = str(self.start)
        l[4] = str(self.end)
        l[6] = self.strand
        return l

    
    def __str__(self):
        if self.type == "header":
            return " ".join(self.get_line()) +"\n"
        else:
            return "\t".join(self.get_line()) + "\n"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return self.type == other
    
    def get_name(self):
        return self.description["Name"]


def parse_gff(path):
    gff = dict()

    with open(path) as ingff:
        for line in ingff:
            l = line.strip().split("\t")
            
            entry = gff_entry(l)

            if entry == "region" or entry == "header":
                if entry.type in gff:
                    gff[entry.type + ".1"] = entry
                    continue
                gff[entry.type] = entry
                continue

            gff[(entry.type, entry.get_name())] = entry
    return gff

def parse_fasta(path):
    fasta_dict = defaultdict(str)
    with open(path, "r") as infasta:
        for line in infasta:
            if line.startswith(">"):
                ind = line.strip()
                continue
            fasta_dict[ind]+= line.strip()
    return fasta_dict

def write_fasta(fastadict, path):
    with open(path, "w") as outf:
        outf.write("".join([f"{x}\n{fastadict[x]}\n" for x in fastadict]))

def write_gff(gff, path):
    outgff = open(path, "w")
    for i in gff:
        outgff.write(str(gff[i]))
    outgff.close()
        
def shift_pos(fasta, pos):
    print(pos)
    return {x: fasta[x][pos:] + fasta[x][:pos] for x in fasta}

def shift_gff_coord(gff, pos, seq_length):
    for g in gff:

        if gff[g] == "header" or gff[g] =="region":
                continue
        if gff[g].start < pos < gff[g].end:
            print(f"Shifting coordinate would break gene {g}. Aborting")
            exit()

        gff[g].start = gff[g].start - pos 
        gff[g].end = gff[g].end - pos      

        if gff[g].start < 0:
            gff[g].start += seq_length
            gff[g].end += seq_length


def get_gene_pos(gff, gene):
    
    gene_pos = [gff[x].start for x in gff if x != "header" and  gff[x].get_name() == gene]
    if len(gene_pos) == 0:
        print("Missing gene")
        exit()

    return gene_pos[0]

def revcon(sequence):
    return sequence.translate(sequence.maketrans('ACGTacgtRYMKrymkVBHDvbhdNn-', 'TGCAtgcaYRKMyrkmBVDHbvdhNn-'))[::-1]

def reverse_complement_fasta(fasta):
    return {x: revcon(fasta[x]) for x in fasta}

def reverse_complement_gff(gff, genome_length):
    for g in gff:
        if gff[g] == "header" or gff[g] == "region":
            continue
        
        gff[g].strand = "+" if gff[g].strand == "-" else "-"
        s,e = gff[g].start, gff[g].end
    
        gff[g].end = genome_length + 1 - s
        gff[g].start = genome_length + 1 - e

def extract_coding_seq(fasta, gff):
    exons =  [gff[x] for x in gff if gff[x] == "exon"]
    genes = [gff[x] for x in gff if gff[x] == "gene"]
    seqdict = {}
    
    for g in genes:
        child_exons= [x for x in exons if x.get_name() == g.get_name()] 
        seq = ""
        for exon in child_exons:
            start = exon.start
            end = exon.end
            s = exon.strand
            if s == "+":
                seq += fasta[(start -1):end]
            else:
                seq+= revcon(fasta[(start -1):end])
        seqdict[">" + g.get_name()] = seq
    return seqdict

def compare_fasta(fasta1, fasta2):
    if any([x not in fasta2 for x in fasta1]) or any([x not in fasta1 for x in fasta2]):
        print("fasta have different entry")
        return False

    if any([fasta1[x] != fasta2[x] for x in fasta1]):
        return False
    return True

def main():
    fasta_path, gff_path, pos, outprefix = sys.argv[1:5]
    
    print(f"parsing files: {fasta_path} {gff_path}")
    fasta = parse_fasta(fasta_path)
    gff = parse_gff(gff_path)
    genome_length = len(list(fasta.values())[0])
    
    print(f"extracting coding sequences for quality check")
    original_CDS = extract_coding_seq(list(fasta.values())[0], gff)
    
    
    print("checking if require reverse complementing sequence")

    if gff[("gene", "cox1")].strand == "-":
         print("Reverse-complementing sequences and annotations")
         new_fasta = reverse_complement_fasta(fasta)
         reverse_complement_gff(gff, genome_length)

         print("CDS sequence check..")
         new_CDS = extract_coding_seq(list(new_fasta.values())[0], gff)

         if not compare_fasta(new_CDS, original_CDS):
             print("ERROR, CDS modified by script")

    try:
        target_pos = int(pos)
    except:
        target_pos = get_gene_pos(gff, pos) - 1

    print(f"shifting positions by {pos}")
    new_fasta = shift_pos(new_fasta, target_pos)
    shift_gff_coord(gff, target_pos, genome_length)
    
    print(f"extracting new CDS sequence for comparison")
    new_CDS = extract_coding_seq(list(new_fasta.values())[0], gff)
    
    if not compare_fasta(new_CDS, original_CDS):
        print("ERROR, CDS modified by script")
        exit()
                
    write_fasta(new_fasta, outprefix + ".fasta")
    write_gff(gff, outprefix + ".gff")
    


if __name__ == "__main__":
    main()

