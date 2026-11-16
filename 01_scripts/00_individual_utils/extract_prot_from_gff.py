#!/bin/python3

# Modules
import sys, itertools, argparse

# User input
p = argparse.ArgumentParser()
p.add_argument("-f", required = True, help = "Path to fasta file")
p.add_argument("-g", required = True, help = "Path to gff3 file" )
p.add_argument("-taxtable", default = None, required = False, help = "NCBI-entry to SP name correspondance table")
p.add_argument("-blacklist", default = None, required = False, help = "path to a file listing sequences from fasta to ignore")
args = p.parse_args()

if args.blacklist is not None:
    blacklist_path = args.blacklist
    with open(blacklist_path, "r") as inf:
        args.blacklist = [x.strip().replace(" ", "_") for x in inf]
else:
    args.blacklist = []

#Constant
OUTDIR = "15_new_gene_sequences"
gene_names = {
    "ATP6": "ATP6", "atp6": "ATP6",
    "ATP8": "ATP8", "atp8": "ATP8",
    "COI": "COI", "COX1": "COI", "cox1": "COI", "COXI": "COI",
    "COII": "COII", "COX2": "COII", "cox2": "COII",
    "COIII": "COIII", "COX3": "COIII", "cox3": "COIII",
    "CYTB": "CYTB", "cob": "CYTB", "Cytb": "CYTB", "cytb": "CYTB",
    "NADH1": "NADH1","nad1": "NADH1", "ND1": "NADH1",
    "NADH2": "NADH2", "nad2": "NAH2", "NAD2": "NADH2", "ND2": "NADH2",
    "NADH3": "NADH3", "nad3": "NADH3", "ND3": "NADH3",
    "NADH4": "NADH4", "nad4": "NADH4", "ND4": "NADH4",
    "NADH4L": "NADH4L", "nad4l": "NADH4L", "ND4L": "NADH4L",
    "NADH5": "NADH5", "nad5": "NADH5", "ND5": "NADH5",
    "NADH6": "NADH6", "nad6": "NADH6", "ND6": "NADH6"
}

# Functions
def translate_codon(codon):
    Ambiguity_table = {"R": ["A", "G"],
                       "Y": ["C", "T"],
                       "M": ["A", "C"],
                       "K": ["G", "T"],
                       "S": ["G", "C"],
                       "W": ["T", "A"],
                       "B": ["T", "C", "G"],
                       "D": ["A", "T", "G"],
                       "H": ["A", "C", "T"],
                       "V": ["A", "C", "G"],
                       "N": ["A", "T", "C", "G"],
                       "T": "T", "A":"A", "C": "C", "G": "G"}
    
    codontable = {
                    'ATA':'M', 'ATC':'I', 'ATT':'I', 'ATG':'M',
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
    
    if codon in codontable:
        return codontable[codon]

    possible_codons = list(itertools.product(*[Ambiguity_table[c] for c in codon]))
    possible_AA = [codontable["".join(x)] for x in possible_codons]

    if len(set(possible_AA)) == 1:
        return possible_AA[0]

    print(f"Error parsing codon: {codon}; unable to resolve ambuiguity")
    return "X"
        
def translate_seq(seq):
    codons = [seq[x:x+3] for x in range(0, len(seq), 3)]

    if codons[0] in ["ATT", "GTG", "TTG", "ATC"]:
        print(f"Using {codons[0]} alternative start codon")
        codons[0] = "ATG"
    trans_codons = [translate_codon(x) for x in codons if len(x) == 3]
    
    if trans_codons[-1] == "_":
        trans_codons.pop()

    return "".join(trans_codons)

def get_Name(gff_line):
    field = [x for x in gff_line[8].split(";") if x.startswith("Name=")]
    if len(field) == 0:
        return None
    return field[0].split("=")[1]

def conv_acc(acc, tax):
    n = acc[:-2] if acc.endswith(".1") else acc
    if n in tax:
        return tax[n]
    return n

def parse_fasta(path, tax):
    seq_dict = {}

    with open(path, "r") as fastin:
        for line in fastin:
            if line.startswith(">"):
                name = line.strip().split()[0][1:]
                name = conv_acc(name, tax)

                seq_dict[name] = ""
            else:
                seq_dict[name] += line.strip()
    return seq_dict

def parse_gff(path, gene_names, tax):
    genes = {}

    with open(path, "r") as gffin:
        for x, line in enumerate(gffin):
            if line.startswith("#") or line.strip() == "":
                continue
            l = line.strip().split()
            if l[2] != "gene":
                continue
            name = get_Name(l)
            
            if name is None:
                continue
            
            if name not in gene_names:
                print(f"Issue parsing: {l}; {name} ")
                continue
            
            name_clr = gene_names[name]
            
            sp_name = conv_acc(l[0], tax)
            genes[sp_name + "-" + name_clr] = [int(l[3]), int(l[4]), l[6]]

    return genes
    
def parse_taxtable(table):
    convtable = {}
    if table is None:
        return convtable

    with open(table, "r") as tablin:
        for line in tablin:
            print(line)
            l = line.strip().split()
            acc = l[0]
            if acc.endswith(".1"):
                acc = acc[:-2]
            convtable[acc] = l[1]
    return convtable

def parse_gene_seq(coord, fasta, blacklist):
    results_prot = {} 
    for entry in coord:
        start, end, strand = coord[entry]
        sp, gene = entry.split("-")
        if sp in blacklist:
            continue
        seq = fasta[sp][(start - 1):end]
        if strand == "-":
            seq = seq.translate(seq.maketrans('ACGTacgtRYMKrymkVBHDvbhdNn-', 'TGCAtgcaYRKMyrkmBVDHbvdhNn-'))[::-1]
        prot = translate_seq(seq)
        if gene not in results_prot:
            results_prot[gene] = ""
        results_prot[gene] += f">{sp}\n{prot}\n"
    return results_prot


# Main
tax = parse_taxtable(args.taxtable)
gene_coord = parse_gff(args.g, gene_names, tax)
fasta_dict = parse_fasta(args.f, tax)
prot_seq = parse_gene_seq(gene_coord, fasta_dict, args.blacklist)

for gene in prot_seq:
    with open(f"{OUTDIR}/{gene}.prot", "w") as prot:
        prot.write(prot_seq[gene])

