#!/bin/bash
# the reference genome sequence used for read mapping
ref_genome="../data/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full.fna"
#ref_genome="../data/genomes/Klebsiella_pneumoniae_KPNIH1/Klebsiella_pneumoniae_KPNIH1.fna"

reads_1="../read-datasets/KPNIH1/SRR1505904_pass_1.fastq"
reads_2="../read-datasets/KPNIH1/SRR1505904_pass_1.fastq"

out_dir="../read-mapping/kp-kpnih1-back-mutated-full-real/"
file_prefix="kp-back-mutated-full"

# ==== Read mapping ====

# Build index for reference genome
bowtie2-build $ref_genome $out_dir"genome-index"

# Single mapping, best-match for MTB
bowtie2 -x $out_dir"genome-index" -U $reads_1,$reads_2 -S $out_dir$file_prefix-mapping-best-match.sam 2> $out_dir$file_prefix-mapping-best-match-log.txt

# Single mapping, report-all for MTB
bowtie2 -a -x $out_dir"genome-index" -U $reads_1,$reads_2 -S $out_dir$file_prefix-mapping-report-all.sam 2> $out_dir$file_prefix-mapping-report-all-log.txt

echo "\n=== Read mapping DONE! ===\n"

# ==== Multimapping Resolution ====


cd ..
prom="simple-bayesian.py"
remu="bayesian-update-main.py"
out_dir="./read-mapping/kp-kpnih1-back-mutated-full-real/"

python3 $prom > $out_dir$file_prefix-log-prom.txt

echo "\n=== PROM DONE! ===\n"

python3 $remu > $out_dir$file_prefix-log-remu.txt

echo "\n=== REMU DONE! ===\n"

cd bash-scripts

# ==== SAMTools ====

out_dir="../read-mapping/kp-kpnih1-back-mutated-full-real/"

# Creating BAM files for Bowtie2 best-match
samtools view -bS $out_dir$file_prefix-mapping-best-match.sam | samtools sort - -o $out_dir$file_prefix-mapping-best-match-sorted.bam
samtools index $out_dir$file_prefix-mapping-best-match-sorted.bam

# Creating BAM files for Bowtie2 report-all
samtools view -bS $out_dir$file_prefix-mapping-report-all.sam | samtools sort - -o $out_dir$file_prefix-mapping-report-all-sorted.bam
samtools index $out_dir$file_prefix-mapping-report-all-sorted.bam

# Creating BAM files for PROM
samtools view -bS $out_dir$file_prefix-mapping-prom.sam | samtools sort - -o $out_dir$file_prefix-mapping-prom-sorted.bam
samtools index $out_dir$file_prefix-mapping-prom-sorted.bam

# Creating BAM files for REMU
samtools view -bS $out_dir$file_prefix-mapping-remu.sam | samtools sort - -o $out_dir$file_prefix-mapping-remu-sorted.bam
samtools index $out_dir$file_prefix-mapping-remu-sorted.bam

sh variant_calling.sh

echo "\n=== SAMTools DONE! ===\n"
