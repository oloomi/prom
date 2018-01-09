#!/bin/bash

# the genome from where synthetic reads are generated
reads_genome="./genome-mutated/mutated-genome.fna"

# the reference genome sequence used for read mapping
ref_genome="../../genome-ref/ref-genome.fna"

# ==== Generating synthetic reads ====

reads_out_dir="./reads/"
reads_file_prefix="reads"

art="art_illumina"
#Simulation of single-end reads of 150 bp with coverage 25; with max number of indels = 0 (-k 0)
$art -ss HS25 -sam -i $reads_genome -l 150 -f 25 -k 0 -rs 12345 -o $reads_out_dir$reads_file_prefix

# Creating BAM files for ArtIllumina benchmark
#samtools view -bS $reads_out_dir$reads_file_prefix.sam -o $reads_out_dir$reads_file_prefix.bam
#samtools sort $reads_out_dir$reads_file_prefix.bam -o $reads_out_dir$reads_file_prefix-sorted.bam
#samtools index $reads_out_dir$reads_file_prefix-sorted.bam

echo "\n==== Generating synthetic reads DONE! ====\n"

# ==== Bowtie2 Read mapping ====

out_dir="./mappings/bowtie/"
file_prefix="bowtie"

# Build index for reference genome
bowtie2-build $ref_genome $out_dir"genome-index"

# Single mapping, best-match
bowtie2 -x $out_dir"genome-index" -U $reads_out_dir$reads_file_prefix.fq -S $out_dir$file_prefix-mapping-best-match.sam \
2> $out_dir$file_prefix-mapping-best-match-log.txt

# Single mapping, report-all
bowtie2 -a -x $out_dir"genome-index" -U $reads_out_dir$reads_file_prefix.fq -S $out_dir$file_prefix-mapping-report-all.sam \
2> $out_dir$file_prefix-mapping-report-all-log.txt

echo "\n=== Bowtie2 read mapping DONE! ===\n"

# ==== SAMTools on Bowtie2 ====

# Creating BAM files for best-match
samtools view -bS $out_dir$file_prefix-mapping-best-match.sam -o $out_dir$file_prefix-mapping-best-match.bam
samtools sort $out_dir$file_prefix-mapping-best-match.bam -o $out_dir$file_prefix-mapping-best-match-sorted.bam
samtools index $out_dir$file_prefix-mapping-best-match-sorted.bam

# Creating BAM files for report-all
samtools view -bS $out_dir$file_prefix-mapping-report-all.sam -o $out_dir$file_prefix-mapping-report-all.bam
samtools sort $out_dir$file_prefix-mapping-report-all.bam -o $out_dir$file_prefix-mapping-report-all-sorted.bam
samtools index $out_dir$file_prefix-mapping-report-all-sorted.bam

echo "\n=== Bowtie 2 SAMTools DONE! ===\n"

# ==== BWA Read mapping ====

out_dir="./mappings/bwa/"
file_prefix="bwa"
cd $out_dir
out_dir="./"
ref_genome="../../../../genome-ref/ref-genome.fna"
reads_out_dir="../../reads/"

# Build index for reference genome
bwa index $ref_genome genome-index

# Single mapping, best-match
bwa mem genome-index $reads_out_dir$reads_file_prefix.fq > $out_dir$file_prefix-mapping-best-match.sam \
2> $out_dir$file_prefix-mapping-best-match-log.txt

# Single mapping, report-all
bwa mem -a genome-index $reads_out_dir$reads_file_prefix.fq > $out_dir$file_prefix-mapping-report-all.sam \
2> $out_dir$file_prefix-mapping-report-all-log.txt

echo "\n=== BWA read mapping DONE! ===\n"

# ==== SAMTools on BWA ====

# Creating BAM files for best-match
samtools view -bS $out_dir$file_prefix-mapping-best-match.sam -o $out_dir$file_prefix-mapping-best-match.bam
samtools sort $out_dir$file_prefix-mapping-best-match.bam -o $out_dir$file_prefix-mapping-best-match-sorted.bam
samtools index $out_dir$file_prefix-mapping-best-match-sorted.bam

# Creating BAM files for report-all
samtools view -bS $out_dir$file_prefix-mapping-report-all.sam -o $out_dir$file_prefix-mapping-report-all.bam
samtools sort $out_dir$file_prefix-mapping-report-all.bam -o $out_dir$file_prefix-mapping-report-all-sorted.bam
samtools index $out_dir$file_prefix-mapping-report-all-sorted.bam

cd ../..
echo "\n=== BWA SAMTools DONE! ===\n"


