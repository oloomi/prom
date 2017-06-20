#!/bin/bash

#art="E:\BioinfTools\art\art_illumina"
export PATH=$PATH:/home/mohammad/Applications/art_bin_MountRainier/
art="art_illumina"

# the genome from where synthetic reads are generated
#reads_genome="../data/genomes/toy-genome-mutated.fna"
#reads_genome="../data/genomes/ot-whole-genome-mutated-70-140.fna"
#reads_genome="../data/genomes/mtb-whole-genome-mutated-70-140.fna"
reads_genome="../data/genomes/mtb-whole-genome-mutated-100-140-half.fna"
#reads_genome="../data/genomes/mtb-whole-genome-mutated.fna"
#reads_genome="../data/genomes/mtb-genome-extract.fna"
#reads_genome="../data/genomes/mtb-genome-extract-mutated-long-repeats.fna"

# the reference genome sequence used for read mapping
#ref_genome="../data/genomes/toy-genome.fna"
#ref_genome="../data/genomes/Orientia_tsutsugamushi_Ikeda_uid58869/NC_010793.fna"
ref_genome="../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna"
#ref_genome="../data/genomes/mtb-genome-extract.fna"


#out_dir="../read-mapping/toy-genome-mutated/"
#file_prefix="toy-wg-mutated-se"
#out_dir="../read-mapping/ot-whole-genome-mutated-70-140/"
#file_prefix="ot-wg-mutated-se"
out_dir="../read-mapping/mtb-whole-genome-mutated-100-140/"
file_prefix="mtb-wg-mutated-se"
#out_dir="../read-mapping/mtb-whole-genome-mutated-70-140/"
#out_dir="../read-mapping/mtb-whole-genome-mutated/"
#file_prefix="mtb-wg-mutated-se"
#out_dir="../read-mapping/mtb-mutated-long-repeats/"
#file_prefix="mtb-mutated-se"
#out_dir="../read-mapping/mtb-normal/"
#file_prefix="mtb-normal-se"

# ==== Generating synthetic reads ====

#Simulation of single-end reads of 150 bp with coverage 10; with max number of indels = 0 (-k 0)
#$art -ss HS25 -sam -i $reads_genome -l 150 -f 10 -k 0 -rs 12345 -o $out_dir$file_prefix
$art -ss HS25 -sam -i $reads_genome -l 150 -f 25 -k 0 -rs 12345 -o $out_dir$file_prefix

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
samtools sort $out_dir$file_prefix.bam -o $out_dir$file_prefix-sorted.bam
samtools index $out_dir$file_prefix-sorted.bam

# Creating BAM files for MTB best-match
samtools view -bS $out_dir$file_prefix-mapping-best-match.sam -o $out_dir$file_prefix-mapping-best-match.bam
samtools sort $out_dir$file_prefix-mapping-best-match.bam -o $out_dir$file_prefix-mapping-best-match-sorted.bam
samtools index $out_dir$file_prefix-mapping-best-match-sorted.bam

# Creating BAM files for MTB report-all
samtools view -bS $out_dir$file_prefix-mapping-report-all.sam -o $out_dir$file_prefix-mapping-report-all.bam
samtools sort $out_dir$file_prefix-mapping-report-all.bam -o $out_dir$file_prefix-mapping-report-all-sorted.bam
samtools index $out_dir$file_prefix-mapping-report-all-sorted.bam


echo "\n=== SAMTools DONE! ===\n"
