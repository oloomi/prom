#!/bin/bash
# Running other multi-mapping resolution methods

reference="../../genome-ref/ref-genome.fna"

# ----------- Bowtie2 ----------
alignments="./mappings/bowtie/bowtie-mapping-report-all-sorted"
alignments_sam="./mappings/bowtie/bowtie-mapping-report-all"
outfile="./mappings/bowtie/bowtie"

# MMR method
samtools sort -n $alignments.bam -o $alignments-id-sorted.bam

/usr/bin/time -v mmr -o $outfile-mmr.bam -F 3 -b -R 150 $alignments-id-sorted.bam

sam_file=$outfile-mmr
samtools sort $sam_file.bam -o $sam_file-sorted.bam
samtools index $sam_file-sorted.bam

echo "\n=== Bowtie + MMR Multimapping Resolution Done ===\n"

# REMU method
/usr/bin/time -v remu.py -g $reference -i $alignments_sam.sam -o $outfile-remu.sam -r 10

sam_file=$outfile-remu
samtools view -bS $sam_file.sam -o $sam_file.bam
samtools sort $sam_file.bam -o $sam_file-sorted.bam
samtools index $sam_file-sorted.bam

echo "\n=== Bowtie + REMU Multimapping Resolution Done ===\n"

# ----------- BWA ----------
alignments="./mappings/bwa/bwa-mapping-report-all-sorted"
alignments_sam="./mappings/bwa/bwa-mapping-report-all"
outfile="./mappings/bwa/bwa"

# MMR method
samtools sort -n $alignments.bam -o $alignments-id-sorted.bam

/usr/bin/time -v mmr -o $outfile-mmr.bam -F 3 -b -R 150 $alignments-id-sorted.bam

sam_file=$outfile-mmr
samtools sort $sam_file.bam -o $sam_file-sorted.bam
samtools index $sam_file-sorted.bam

echo "\n=== BWA + MMR Multimapping Resolution Done ===\n"

# REMU method
/usr/bin/time -v remu.py -g $reference -i $alignments_sam.sam -o $outfile-remu.sam -r 10

sam_file=$outfile-remu
samtools view -bS $sam_file.sam -o $sam_file.bam
samtools sort $sam_file.bam -o $sam_file-sorted.bam
samtools index $sam_file-sorted.bam

echo "\n=== BWA + REMU Multimapping Resolution Done ===\n"
