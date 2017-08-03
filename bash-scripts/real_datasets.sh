#!/bin/bash

# the reference genome sequence used for read mapping
ref_genome="../data/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full.fna"
#ref_genome="../data/genomes/Klebsiella_pneumoniae_KPNIH1/Klebsiella_pneumoniae_KPNIH1.fna"

reads_1="../read-datasets/KPNIH1/SRR1505904_pass_1.fastq"
reads_2="../read-datasets/KPNIH1/SRR1505904_pass_1.fastq"

# Remember to revert back
out_dir="../read-mapping/kp-kpnih1-back-mutated-full-real/"
file_prefix="kp-back-mutated-full"





# Creating BAM files for PROM
samtools view -bS $out_dir$file_prefix-mapping-prom.sam | samtools sort - -o $out_dir$file_prefix-mapping-prom-sorted.bam
samtools index $out_dir$file_prefix-mapping-prom-sorted.bam

sh variant_calling.sh

echo "\n=== SAMTools DONE! ===\n"
