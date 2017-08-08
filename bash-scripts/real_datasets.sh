#!/bin/bash

ref_genome="../data/genomes/MTB-H37Rv-back-mutated-full.fna"

reads_1="../read-datasets/MTB-H37Rv/fastq/SRR2818101_pass_1_trim_250.fastq"

out_dir="../read-mapping/mtb-h37rv-back-mutated/"
file_prefix="mtb-h37rv-back-mutated"

# ==== Multimapping Resolution ====

cd ..
prom="simple-bayesian.py"
remu="bayesian-update-main.py"
out_dir="./read-mapping/mtb-h37rv-back-mutated/"
#out_dir="./read-mapping/kp-kpnih1-back-mutated-full-real/"

#python3 $prom > $out_dir$file_prefix-log-prom.txt

echo "\n=== PROM DONE! ===\n"

#python3 $remu > $out_dir$file_prefix-log-remu.txt

echo "\n=== REMU DONE! ===\n"

cd bash-scripts

# ==== SAMTools ====

out_dir="../read-mapping/mtb-h37rv-back-mutated/"
#out_dir="../read-mapping/kp-kpnih1-back-mutated-full-real/"

# Creating BAM files for Bowtie2 best-match
samtools view -bS $out_dir$file_prefix-mapping-best-match.sam | samtools sort - -o $out_dir$file_prefix-mapping-best-match-sorted.bam
samtools index $out_dir$file_prefix-mapping-best-match-sorted.bam

# Creating BAM files for Bowtie2 report-all
samtools view -bS $out_dir$file_prefix-mapping-report-all.sam | samtools sort - -o $out_dir$file_prefix-mapping-report-all-sorted.bam
samtools index $out_dir$file_prefix-mapping-report-all-sorted.bam

# Creating BAM files for PROM
#samtools view -bS $out_dir$file_prefix-mapping-prom.sam | samtools sort - -o $out_dir$file_prefix-mapping-prom-sorted.bam
#samtools index $out_dir$file_prefix-mapping-prom-sorted.bam

# Creating BAM files for REMU
#samtools view -bS $out_dir$file_prefix-mapping-remu.sam | samtools sort - -o $out_dir$file_prefix-mapping-remu-sorted.bam
#samtools index $out_dir$file_prefix-mapping-remu-sorted.bam

sh variant_calling.sh

echo "\n=== SAMTools DONE! ===\n"
