# Arctic_Leafhopper_mitochondria_paper

This repos contain all the script necassary to reproduce the results presented in the paper:
### Mitochondrial genomes of leafhoppers and a parasitoid fly expand genomic resources for the Eastern Canadian Arctic
### Sylvestre F., Plante N., Molligan J., Grégoire Taillefer A., Pérez-lopez E. 
### LINK Avaialable upon publication


## Requirements:
- python 3.12
- fastp 0.15.0
- GetOrganelle v 1.7.7.1  (https://github.com/kinggerm/getorganelle)
- Bandage v0.8.1
- tRNAsca-SE v2.0
- MAFFT v7.525
- bwa v0.7.18
- mosdepth v0.3.12
- bcftools v1.21
- https://molevol-ibe.csic.es/Gblocks.html

## Initialisation:
Sample reads can be find on NCBI in Bioproject PRJNA1450552.
We included:
AL2: Diplocolenus abdoninalis (verdanus)
AL3: Hebecephalus algidus
AL5: Macustus grisescens
AL6: Psammotettix alienus
AL7: Limotettix myralis
AL8: limotettix plutonius
AL9: Stroggylocephalus mixtus

Download them in the 02_raw_files folder if you want to run the pipeline

## Mitochondria assembly
- Script 1-2 clean the reads and provide a first assembly of the mitochondria usign GetOrganelle
- These assembly are fragmented due to the presence of strong polymorphism
- Genome graph were curated using Bandage, a manual step that couldn't be reproduced here
- After that, script 3 was use to obtain new mitochondria assembly, which were extracted and renamed
- For species with multiple mitochondria assembled, each was ID using the BOLD portal (https://id.boldsystems.org/), renamed and reID accordingly 
- After this step, put each assembled fasta in 05_sequences_and_annotations
- Run mitos2 on galaxy (https://microbiology.usegalaxy.org/?tool_id=mitos2) using invertabrate code and the RefSeq89 Metazoa db and download resulting gff3 in 05_sequences_and_annotations

## Mitochondria cleaning
- Scripts 4-8 Reorient mitochondria, re-annotate tRNA to complete their annotation and allow manual cleaning of start and stop codon before providing an NCBI-ready gff file
- Scripts 9-11 realign reads back to the associated mitochondria, before estimating coverage and SNP density for visualisation in proksee (https://proksee.ca/)
- Scripts 12-12 extract and align protein-coding sequences from generated mitochondria and NCBI-available leafhopper mitochondria (available in 99_paper_specific_files)

# Phylogenetic analysis
- To reproduce phylogenetic analysis, run G-block on concatenated_alignement.fasta available in 18_aligned_sequences
- then run IQTREE as:

iqtree2 -s 18_aligned_sequences/concatenated_alignement.fasta-gb.fa -B 1000 -MFP








