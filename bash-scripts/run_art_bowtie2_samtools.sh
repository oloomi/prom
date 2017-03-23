#!/bin/bash

#art="E:\BioinfTools\art\art_illumina"
export PATH=$PATH:/home/mohammad/Applications/art_bin_MountRainier/
art="art_illumina"

# the genome from where synthetic reads are generated
#reads_genome="../data/genomes/mtb-genome-extract.fna"
reads_genome="../data/genomes/mtb-genome-extract-mutated.fna"

# the reference genome sequence used for read mapping
ref_genome="../data/genomes/mtb-genome-extract.fna"

#out_dir="../read-mapping/mtb-normal/"
#file_prefix="mtb-normal-se"
out_dir="../read-mapping/mtb-mutated/"
file_prefix="mtb-mutated-se"

# ==== Generating synthetic reads ====

#Simulation of single-end reads of 150 bp with coverage 10; with max number of indels = 0 (-k 0)
$art -ss HS25 -sam -i $reads_genome -l 150 -f 10 -k 0 -rs 12345 -o $out_dir$file_prefix

echo "\n==== Generating synthetic reads DONE! ====\n"

# ==== Read mapping ====

# Build index for MTB reference genome
bowtie2-build $ref_genome $out_dir"genome-index"

# Single mapping, best-match for MTB
bowtie2 -x $out_dir"genome-index" -U $out_dir$file_prefix.fq -S $out_dir$file_prefix-mapping-best-match.sam 2> $out_dir$file_prefix-mapping-best-match-log.txt

# Single mapping, report-all for MTB
bowtie2 -a -x $out_dir"genome-index" -U $out_dir$file_prefix.fq -S $out_dir$file_prefix-mapping-report-all.sam 2> $out_dir$file_prefix-mapping-report-all-log.txt

echo "\n=== Read mapping DONE! ===\n"

# ==== SAMTools ====

# Creating BAM files for MTB ArtIllumina benchmark
samtools view -bS $out_dir$file_prefix.sam -o $out_dir$file_prefix.bam
samtools sort $out_dir$file_prefix.bam $out_dir$file_prefix-sorted
samtools index $out_dir$file_prefix-sorted.bam

# Creating BAM files for MTB best-match
samtools view -bS $out_dir$file_prefix-mapping-best-match.sam -o $out_dir$file_prefix-mapping-best-match.bam
samtools sort $out_dir$file_prefix-mapping-best-match.bam $out_dir$file_prefix-mapping-best-match-sorted
samtools index $out_dir$file_prefix-mapping-best-match-sorted.bam

# Creating BAM files for MTB report-all
samtools view -bS $out_dir$file_prefix-mapping-report-all.sam -o $out_dir$file_prefix-mapping-report-all.bam
samtools sort $out_dir$file_prefix-mapping-report-all.bam $out_dir$file_prefix-mapping-report-all-sorted
samtools index $out_dir$file_prefix-mapping-report-all-sorted.bam


echo "\n=== SAMTools DONE! ===\n"