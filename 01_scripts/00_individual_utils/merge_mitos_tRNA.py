#!/bin/python

import sys
from collections import defaultdict
from pathlib import Path

codontable = {'ATA':'M', 'ATC':'I', 'ATT':'I', 'ATG':'M',
               'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
               'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
               'AGC':'S', 'AGT':'S', 'AGA':'S', 'AGG':'S',
               'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
               'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
               'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
               'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
               'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
               'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
               'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
               'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
               'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
               'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
               'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
               'TGC':'C', 'TGT':'C', 'TGA':'W', 'TGG':'W',
               }


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
            self.description = {x.split("=")[0]:x.split("=")[1] for x in line[8].strip(";").split(";")}

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
        return self.description["Name"] if "Name" in self.description else self.description["Parent"]
    
    def get_anticodon(self):
        if self.type != "tRNA":
            return None

        if "anticodon" in self.description:
            return self.description["anticodon"]
        
        name = self.get_name().split("(")[1][:-1].upper()
        return name


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

def write_gff(gff, path):
    outgff = open(path, "w")
    for i in gff:
        outgff.write(str(gff[i]))
    outgff.close()
        
def check_trna_in_mitos(entry, gff_dict):
    found = False
    for en in gff_dict:
        if gff_dict[en] != "tRNA":
            continue

        if abs(gff_dict[en].start - entry.start) <= 5 and abs(gff_dict[en].end - entry.end) < 5:
            if gff_dict[en].get_anticodon() == entry.get_anticodon():
                found = True

    return found
            
def conv_trna_to_mitos(trna_entry, second):
    inc = "2" if second else ""
    codon = trna_entry.get_anticodon()
    aa = get_aa_anticodon(codon)
    new_description = {"ID": f"transcript_trn{aa}{inc}({codon.lower()})",
                              "Name" : f"trn{aa}{inc}({codon.lower()})",
                              "Parent" : f"gene_trn{aa}{inc}({codon.lower()})",
                              "gene_id" : f"trn{aa}{inc}({codon.lower()})"}
    line = trna_entry.line.copy()
    line[8] = ";".join([f"{x}={new_description[x]}" for x in new_description])
    return gff_entry(line)

def get_aa_anticodon(anticodon):
    codon = anticodon.translate(anticodon.maketrans('ACGTacgtRYMKrymkVBHDvbhdNn-', 'TGCAtgcaYRKMyrkmBVDHbvdhNn-'))[::-1]
    return codontable[codon]

def create_exon_line(trna_entry, second):
    inc = "2" if second else ""
    codon = trna_entry.get_anticodon()
    aa = get_aa_anticodon(codon)
    line = trna_entry.line.copy()
    line[1] = "tRNAscan-SE"
    line[2] = "exon"
    line[8] = f"Parent=transcript_trn{aa}{inc};Name=trn{aa}{inc}"
    return gff_entry(line)
    

def create_transcript_line(trna_entry, second):
    inc = "2" if second else ""
    codon = trna_entry.get_anticodon()
    aa = get_aa_anticodon(codon)
    line = trna_entry.line.copy()
    line[1] = "tRNAscan-SE"
    line[2] = "ncRNA_gene"
    line[8] = f"ID=gene_trn{aa}{inc};Name=trn{aa}{inc};gene_id=trn{aa}{inc}"
    return gff_entry(line)


def main():
    mitos_gffp, trnaseq_gffp = sys.argv[1:3]
    outpath = "08_merged_annotations"
    print(f"parsing files: {mitos_gffp}")
    mitos_gff = parse_gff(mitos_gffp)
    trna_gff = parse_gff(trnaseq_gffp)
    dup = False
    seq_to_add = [trna_gff[x].get_name() for x in trna_gff if trna_gff[x] == "tRNA" and not check_trna_in_mitos(trna_gff[x], mitos_gff)]
    
    for i in seq_to_add:
        for x in trna_gff:
            dup = False
            if trna_gff[x].get_name() != i:
                continue

            if trna_gff[x] == "tRNA":
                if ("ncRNA_gene", f"trn{get_aa_anticodon(trna_gff[x].get_anticodon())}") in mitos_gff:
                    dup = True
                gene = create_transcript_line(trna_gff[x], dup)
                new_trna = conv_trna_to_mitos(trna_gff[x], dup)
                new_exon = create_exon_line(trna_gff[x], dup)


                mitos_gff[("ncRNA_gene", gene.get_name())] = gene
                mitos_gff[("tRNA", new_trna.get_name())] = new_trna
                mitos_gff[("exon", new_exon.get_name())] = new_exon
                
    

    write_gff(mitos_gff, outpath + f"/{Path(mitos_gffp).name}")


if __name__ == "__main__":
    main()
