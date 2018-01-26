#!/bin/bash

make_exp_dir() {
  mkdir genome-mutated mappings reads results variants
  mkdir mappings/bowtie mappings/bwa
}

create_dir() {
  mkdir $1
  cd $1
  mkdir genome-ref real-data simulated-data
  mkdir genome-ref/repeats simulated-data/supermax-100-140 simulated-data/middle-supermax
  cd ./simulated-data/supermax-100-140
  make_exp_dir
  cd ../middle-supermax
  make_exp_dir
  cd ../../read-data
  mkdir back-mutate
  cd ./back-mutate
  make_exp_dir
  cd ../..
}

get_genome() {
  wget -O - $1 | gunzip -c > ./genome-ref/ref-genome.fna
}

find_repeats() {
  cd ./genome-ref/repeats
  mkvtree -db ../ref-genome.fna -dna -v -allout -pl
  vmatch -supermax -l 150 -h 1 ref-genome.fna > supermax-repeats.txt
  vmatch -l 150 -d -p ref-genome.fna > all-repeats.txt
  vmatch -tandem -l 150 ref-genome.fna > tandem-repeats.txt
  cd ../..
}

get_reads() {
  fastq-dump --outdir ./reads/fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-files --clip $1
  gunzip ./reads/fastq/$1_pass_1.fastq.gz
  fastx_trimmer -l 150 -m 150 -Q33 -i ./reads/fastq/$1_pass_1.fastq -o ./reads/reads.fq
}

prepare_data() {
  create_dir $1
  get_genome $2
  find_repeats
}

# Mycobacterium Tuberculosis H37rv
#prepare_data mtb ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz

#fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip SRR2818101
#gunzip fastq/SRR2818101_pass_1.fastq.gz
#fastx_trimmer -l 150 -m 150 -Q33 -i ./fastq/SRR2818101_pass_1.fastq -o ./reads.fq

# Arabidopsis thaliana
#prepare_data athalina ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_genomic.fna.gz

# Fungus (Neurospora crassa), 101 repeats, 40.7 Mb
#prepare_data fungus ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.fna.gz

# Yeast (Saccharomyces cerevisiae), 474 repeats, 12.1 Mb, 16 chromosomes
#prepare_data yeast ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
#real data
#fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip ERR1938683
make_exp_dir
get_reads ERR1938683


# Human chromosome 19 GRCh38.p7 assembly
#prepare_data human-chr19 http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr19.fa.gz

# Ecoli (Escherichia coli), 17 repeats, 5.1 Mb
#prepare_data ecoli ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz


#├── genome-ref
#│   └── repeats
#├── real-data
#└── simulated-data
#    └── supermax-100-140
#        ├── genome-mutated
#        ├── mappings
#        │   ├── bowtie
#        │   └── bwa
#        ├── reads
#        ├── results
#        └── variants
