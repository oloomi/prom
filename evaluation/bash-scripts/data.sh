#!/bin/bash

# MTB
mkdir genome-ref
cd genome-ref
wget -O - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz | gunzip -c > ref-genome.fna

mkdir repeats
cd repeats
mkvtree -db ../ref-genome.fna -dna -v -allout -pl
vmatch -supermax -l 150 ref-genome.fna > supermax-repeats.txt

#vmatch -d -p -l 150 -h 1 NC_000962.fna > mtb-repeats-hamming-1.txt


