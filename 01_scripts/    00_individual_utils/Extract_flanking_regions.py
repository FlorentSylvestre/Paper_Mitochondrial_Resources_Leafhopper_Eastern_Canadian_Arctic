#!/bin/python

"""
Script that extract gene flanking region for inspection of start and stop codon in Geneious

Usage: python3 Extract_gene_flanking_regions <fasta_path>, <gff_path>,  <windows>, <outprefix>
fasta_path: path to a file containing a list of path to fastas, one per line
gff_path: path to a file containing a list of path to gff, one per line
windows: size of gene flanking region to include
outprefix: output directory name
"""


import sys
from collections import defaultdict
from pathlib import Path
class gff_entry:

    def __init__(self, line: None):
        self.line = line

        if line[0].startswith("#"):
            self.type = "header"
            self.description = {"Name":"header"}
            self.species = "NA"
            
        else:
            self.species = line[0]
            self.type = line[2]
            self.start = int(line[3])
            self.end = int(line[4])
            self.strand = line[6]
            self.description = {x.split("=")[0]:x.split("=")[1] for x in line[8].split(";")}

    def get_line(self):
        if self == "header":
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
    
    def get_species(self):
        return self.species

def parse_gff(path):
    gff = dict()
    with open(path) as ingff:
        for line in ingff:
            l = line.strip().split("\t")
            
            entry = gff_entry(l)
            gff[(entry.get_species(), entry.type, entry.get_name())] = entry
    return gff

def parse_fasta(path):
    fasta_dict = defaultdict(str)
    with open(path.strip(), "r") as infasta:
        for line in infasta:
            if line.startswith(">"):
                ind = line.strip()
                continue
            fasta_dict[ind] += line.strip()
    return fasta_dict

def write_fasta(fastadict, path):
    for gene in fastadict:
        with open(path  + f"/{gene}.fasta", "w") as outf:
            outf.write("".join([f"{x}\n{fastadict[gene][x]}\n" for x in fastadict[gene]]))

def write_gff(gff, path):
    for gene in gff:
        outgff = open(path + f"/{gene}.gff", "w")
        for i in gff[gene]:
            outgff.write(str(gff[gene][i]))
        outgff.close()
        
def extract_seq(fasta, gff, windows):

    genes = [x for x in gff if gff[x] == "gene"]
    seqdict = {}
    
    for species, n, g in genes:
        start, end = max(0, gff[species, n, g].start - windows -1), gff[species,n,g].end + windows
        seq = fasta[f">{species}"][start:end]
        if g not in seqdict:
            seqdict[g] = {}
        seqdict[g].update({f">{species}": seq})
    
    return seqdict

def extract_gff(gff, windows):
    genes = [x for x in gff if gff[x] == "gene"]
    gffdict = {}

    for species, n, g in genes:
        start, end = max(1, (gff[species, n, g].start - windows)), gff[species,n,g].end + windows
        sp_gff = {(sp, n, g): gff[(sp, n, g)] for sp, n, g in gff if sp == species}
        entry_kept = [gff_entry(gff[x].line) for x in sp_gff if (start < sp_gff[x].start and end > sp_gff[x].end) or (start > sp_gff[x].start and start < sp_gff[x].end) or (end > sp_gff[x].start and end < sp_gff[x].end)]
        

        for x in entry_kept:
            x.start = max(0, (x.start - start)) + 1
            x.end = min((x.end - start), (end - start)) + 1

        if g not in gffdict:
            gffdict[g] = {}
        gffdict[g].update({(x.get_species(), x.type, x.get_name()): x for x  in entry_kept})
    
    return gffdict

def main():
    fasta_path, gff_path,  windows, outprefix  = sys.argv[1:5]

    print(f"parsing files: {fasta_path}")
    with open(fasta_path, "r") as fasta_list:
        fastas = fasta_list.readlines()

    seqdict = {}
    for fa in fastas:
        seqdict.update(parse_fasta(fa))
    
    gff_list = [Path(gff_path) / Path(x.strip()).with_suffix(".gff").name for x in fastas]
    gff_dict = {}
    
    for gf in gff_list:
        gff_dict.update(parse_gff(gf))
    
    seq = extract_seq(seqdict, gff_dict, int(windows))        
    extracted_gff = extract_gff(gff_dict, int(windows)) 
    
    write_fasta(seq, outprefix)
    write_gff(extracted_gff, outprefix)
    


if __name__ == "__main__":
    main()
