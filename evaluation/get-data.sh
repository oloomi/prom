#!/bin/bash

make_exp_dir() {
  mkdir -p $1/genome-mutated $1/mappings/bowtie $1/mappings/bwa $1/reads $1/results $1/variants
}

create_dir() {
  mkdir $1
  cd $1
  mkdir -p genome-ref/repeats simulated-data/begin-supermax simulated-data/middle-supermax real-data/back-mutate
  make_exp_dir simulated-data/begin-supermax
  make_exp_dir simulated-data/middle-supermax
  make_exp_dir real-data/back-mutate
}

get_genome() {
  wget -O - $1 | gunzip -c > ./genome-ref/ref-genome.fna
}

find_repeats() {
  cd ./genome-ref/repeats
  mkvtree -db ../ref-genome.fna -dna -v -allout -pl
  vmatch -supermax -l $1 -h 1 ref-genome.fna > supermax-repeats.txt
  vmatch -l $1 -d -p ref-genome.fna > all-repeats.txt
  vmatch -tandem -l $1 ref-genome.fna > tandem-repeats.txt
  cd ../..
}

# run_id read_len dir
get_reads() {
  fastq-dump --outdir $3/reads/fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-files --clip $1
  gunzip $3/reads/fastq/$1_pass_1.fastq.gz
  fastx_trimmer -l $2 -m $2 -Q33 -i $3/reads/fastq/$1_pass_1.fastq -o $3/reads/reads.fq
}

# run_id read_len dir url
get_reads_ena() {
  wget $4 -P $3/reads/fastq
  gunzip -c $3/reads/fastq/$1_1.fastq.gz
  fastx_trimmer -l $2 -m $2 -Q33 -i $3/reads/fastq/$1_pass_1.fastq -o $3/reads/reads.fq
}

prepare_data() {
  create_dir $1
  get_genome $3
  find_repeats $2
  cd ..
}

# Mycobacterium Tuberculosis H37rv
prepare_data mtb 150 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz
get_reads SRR2818101 150 mtb/real-data/back-mutate

# Escherichia coli
prepare_data ecoli 100 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
#get_reads ERR022075 100 yeast/real-data/back-mutate
get_reads_ena ERR022075 100 yeast/real-data/back-mutate ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022075/ERR022075_1.fastq.gz
mv yeast/real-data/back-mutate/reads.fq yeast/real-data/back-mutate/all-reads.fq
seqtk sample -s100 yeast/real-data/back-mutate/all-reads.fq 2500000 > yeast/real-data/back-mutate/reads.fq

# Saccharomyces cerevisiae S288C (baker's yeast)
prepare_data yeast 150 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
#get_reads ERR1938683 150 yeast/real-data/back-mutate
get_reads_ena ERR022075 150 yeast/real-data/back-mutate ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR193/ERR1938683/ERR1938683_1.fastq.gz

# Human chromosome 19 GRCh38.p7 assembly
#prepare_data human-chr19 100 http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr19.fa.gz
