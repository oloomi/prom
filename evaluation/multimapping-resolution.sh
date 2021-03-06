#!/bin/bash
# Running other multi-mapping resolution methods

# sh multimapping-resolution.sh [-s -r] [read_len]

if [ "$1" = "-s" ]; then
  ref_genome="../../genome-ref/ref-genome.fna"
elif [ "$1" = "-r" ]; then
  ref_genome="./genome-mutated/mutated-genome.fna"
else
  echo "Invalid argument in variant-calling.sh!"
fi

# ----------- Bowtie2 ----------
alignments="./mappings/bowtie/bowtie-mapping-report-all-sorted"
alignments_sam="./mappings/bowtie/bowtie-mapping-report-all"
outfile="./mappings/bowtie/bowtie"

# MMR method
samtools sort -n ${alignments}.bam -o ${alignments}-id-sorted.bam

/usr/bin/time -v -o ${outfile}-mmr-time-log.txt mmr -o ${outfile}-mmr.bam -F 3 -b -R $2 ${alignments}-id-sorted.bam | tee ${outfile}-mmr-log.txt

sam_file=${outfile}-mmr
samtools sort ${sam_file}.bam -o ${sam_file}-sorted.bam
samtools index ${sam_file}-sorted.bam

echo "\n=== Bowtie + MMR multi-mapping resolution completed! ===\n"

# PROM method
/usr/bin/time -v -o ${outfile}-prom-time-log.txt prom.py -g ${ref_genome} -i ${alignments_sam}.sam -o ${outfile}-prom.sam -r 10 | tee ${outfile}-prom-log.txt

sam_file=${outfile}-prom
samtools view -bS ${sam_file}.sam -o ${sam_file}.bam
samtools sort ${sam_file}.bam -o ${sam_file}-sorted.bam
samtools index ${sam_file}-sorted.bam

echo "\n=== Bowtie + PROM multi-mapping resolution completed! ===\n"

# ----------- BWA ----------
alignments="./mappings/bwa/bwa-mapping-report-all-sorted"
alignments_sam="./mappings/bwa/bwa-mapping-report-all"
outfile="./mappings/bwa/bwa"

# MMR method
samtools sort -n ${alignments}.bam -o ${alignments}-id-sorted.bam

/usr/bin/time -v -o ${outfile}-mmr-time-log.txt mmr -o ${outfile}-mmr.bam -F 3 -b -R $2 ${alignments}-id-sorted.bam | tee ${outfile}-mmr-log.txt

sam_file=${outfile}-mmr
samtools sort ${sam_file}.bam -o ${sam_file}-sorted.bam
samtools index ${sam_file}-sorted.bam

echo "\n=== BWA + MMR multi-mapping resolution completed! ===\n"

# PROM method
/usr/bin/time -v -o ${outfile}-prom-time-log.txt prom.py -g ${ref_genome} -i ${alignments_sam}.sam -o ${outfile}-prom.sam -r 10 | tee ${outfile}-prom-log.txt

sam_file=${outfile}-prom
samtools view -bS ${sam_file}.sam -o ${sam_file}.bam
samtools sort ${sam_file}.bam -o ${sam_file}-sorted.bam
samtools index ${sam_file}-sorted.bam

echo "\n=== BWA + PROM multi-mapping resolution completed! ===\n"
