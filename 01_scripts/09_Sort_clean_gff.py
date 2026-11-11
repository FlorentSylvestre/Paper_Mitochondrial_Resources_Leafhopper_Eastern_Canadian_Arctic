#!/bin/python

import sys
from collections import defaultdict
from pathlib import Path

product_table = {
    "cob": "Cytochrome b", 
    "cox1": "Cytochrome c oxidase, subunit 1",
    "cox2": "Cytochrome c oxidase, subunit 2",
    "cox3": "Cytochrome c oxidase, subunit 3",
    "nad1": "NADH dehydrogenase, subunit 1",
    "nad2": "NADH dehydrogenase, subunit 2",
    "nad3": "NADH dehydrogenase, subunit 3",
    "nad4": "NADH dehydrogenase, subunit 4",
    "nad4l": "NADH dehydrogenase, subunit 4l",
    "nad5": "NADH dehydrogenase, subunit 5",
    "nad6": "NADH dehydrogenase, subunit 6",
    "atp6": "ATP synthase, Fo subunit 6",
    "atp8": "ATP synthase, Fo subunit 8",
    "trnV(tac)": "tRNA-Val",
    "trnL1(tag)": "tRNA-Leu",
    "trnS2(tga)": "tRNA-Ser",
    "trnP(tgg)": "tRNA-Pro",
    "trnT(tgt)": "tRNA-Thr",
    "trnH(gtg)": "tRNA-His",
    "trnF(gaa)": "tRNA-Phe",
    "trnE(ttc)": "tRNA-Glu",
    "trnS1(gct)": "tRNA-Ser",
    "trnS1(tct)": "tRNA-Ser",
    "trnN(gtt)": "tRNA-Asn",
    "trnR(tcg)": "tRNA-Arg",
    "trnA(tgc)": "tRNA-Ala",
    "trnG(tcc)": "tRNA-Gly",
    "trnD(gtc)": "tRNA-Asp",
    "trnK(ctt)": "tRNA-Lys",
    "trnL2(taa)": "tRNA-Leu",
    "trnY(gta)": "tRNA-Tyr",
    "trnC(gca)": "tRNA-Cys",
    "trnW(tca)": "tRNA-Trp",
    "trnM(cat)": "tRNA-Met",
    "trnQ(ttg)": "tRNA-Gln",
    "trnI(gat)": "tRNA-Ile",
    "rrnS": "Ribosomal sub-unit 12S",
    "rrnL": "Ribosomal sub-unit 16S"
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
        if self == "header":
            return self.line 
        l = self.line
        l[2] = str(self.type)
        l[3] = str(self.start)
        l[4] = str(self.end)
        l[6] = self.strand
        l[8] = ";".join([f"{x}={self.description[x]}" for x in self.description])
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
    
    def __lt__(self, other):
        if self.type == "header":
            if other.type != "header":
                return True
            return self.line[0].startswith("##")
        if other.type == "header":
            return False
        if self.type == "region":
            return False
        if other.type == "region":
            return False

        if self.start == other.start:
            if self.type in ["gene", "ncRNA_gene"]:
               return True
            if other.type in ["gene", "ncRNA_gene"]:
               return False
            if self.type in ["mRNA", "rRNA", "tRNA"]:
                return True
            if other.type in ["mRNA", "rRNA", "tRNA"]:
                return False
            if self.type == "exon":
                return False
            if other.type == "exon":
                return True

        if self.start < other.start:
            return True

    def get_name(self):
        return self.description["Name"] if "Name" in self.description else self.description["Parent"]
    

def parse_gff(path):
    gff = list()
    with open(path) as ingff:
        for line in ingff:
            l = line.strip().split("\t")
            gff.append(gff_entry(l))
    return gff

def write_gff(gff, path):
    outgff = open(path, "w")
    for i in gff:
        outgff.write(str(i))
    outgff.close()

def clean_region(gff):
    for y, x in enumerate(gff):
        if x.type != "region":
            continue
        a = str(x.description["Is_circular"]).lower()
        gff[y].description["Is_circular"] = str(x.description["Is_circular"]).lower()
    return gff 

def clean_tRNA_naming(gff):
    for y, x in enumerate(gff):
        if x.type != "tRNA":
            continue
        gff[y-1].description["Name"] = x.description["Name"]
        gff[y-1].description["ID"] = "gene_" + x.description["Name"]
        gff[y-1].description["gene_id"] = x.description["Name"]

        gff[y+1].description["Parent"] = "transcript_" + x.description["Name"]
        gff[y+1].description["Name"] = "transcript_" + x.description["Name"]
    return gff
        
def add_gene_transcript(gff):
    gff_copy = gff.copy()
    shift = 0
    for y, x in enumerate(gff):
        if x.type != "gene":
            continue
        new_line = gff_entry(str(x).split("\t"))
        new_line.type = "mRNA"
        new_line.description["ID"] = "transcript_" + x.description["ID"][5:]
        new_line.description["Parent"] = x.description["ID"]
        new_line.description["Name"] = new_line.description["gene_id"] = x.description["ID"][5:]
        gff_copy.insert( y + shift + 1, new_line)
        shift += 1

    return gff_copy

def add_region_header(gff):
    for y,x in enumerate(gff):
         if x.type != "region":
             continue
         if gff[y-1].line[0].startswith("##"):
             continue
         new_header = f"##sequence-region\t{x.description['ID'].split(":")[0]}\t{x.start}\t{x.end}".split("\t")
         pos = y
    gff.insert(pos, gff_entry(new_header))
    return gff

def remove_OH(gff):
    new_gff = []
    for y, x in enumerate(gff):
        if x.type == "origin_of_replication":
            continue
        
        if x.type == "exon":
            if x.description["Name"].startswith("OH"):
                continue

        new_gff.append(x)

    return new_gff

def add_CDS(gff):
    new_gff = []
    for x, y in enumerate(gff):
        new_gff.append(y)
        if y.type == "exon" and "rrn" not in y.description["Name"] and "trn" not in y.description["Name"]:
            new_entry = gff_entry(str(y).strip().split("\t"))
            new_entry.type = "CDS"
            new_entry.description["Parent"] = "gene-" + y.description["Name"].strip()
            new_gff.append(new_entry)
    return new_gff

def add_gene_descriptor(gff):
    for x in gff:
        if x.type != "gene" and x.type != "ncRNA_gene":
            continue
        x.description["gene"] = x.description["Name"].strip()
    return gff

def add_product(gff):
    for x, y in enumerate(gff):
        if y.type not in ["CDS", "rRNA", "tRNA"]:
            continue
        y.description["product"] = product_table[y.description["Name"].strip()]
    return gff

def clean_tRNA_gene_name(gff):
    for x, y in enumerate(gff):
        if y.type != "ncRNA_gene":
            continue
        y.type = "gene"
    return gff

def main():
    gff_path = sys.argv[1]
    outpath = ".".join(gff_path.split(".")[:-1]) + ".sorted.gff"
    
    print(f"parsing files: {gff_path}")
    gff = parse_gff(gff_path)
    gff.sort()
    write_gff(clean_tRNA_gene_name(
            add_product(
            add_gene_descriptor(
            add_CDS(
            remove_OH(
            add_region_header(
            add_gene_transcript(
            clean_tRNA_naming(
            clean_region(
                gff))))))))),
            outpath)

if __name__ == "__main__":
    main()
