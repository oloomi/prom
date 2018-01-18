#!/bin/bash

create_dir() {
  mkdir genome-ref real-data simulated-data
  mkdir genome-ref/repeats
  mkdir genome-mutated mappings reads results variants
  mkdir mappings/bowtie mappings/bwa
}

get_genome() {
  wget -O - $1 | gunzip -c > ./genome-ref/ref-genome.fna
}

find_repeats() {
  cd ./genome-ref/repeats
  mkvtree -db ../ref-genome.fna -dna -v -allout -pl
  vmatch -supermax -l 150 ref-genome.fna > supermax-repeats.txt
  vmatch -tandem -l 150 -d -p ref-genome.fna > all-repeats.txt
  vmatch -tandem -l 150 ref-genome.fna > tandem-repeats.txt
  cd ../..
}

# Mycobacterium Tuberculosis H37rv
#mkdir mtb
#cd mtb
#create_dir
#get_genome ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz
#find_repeats

# Arabidopsis thaliana
mkdir athaliana
cd athaliana
create_dir
get_genome ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_genomic.fna.gz
find_repeats
