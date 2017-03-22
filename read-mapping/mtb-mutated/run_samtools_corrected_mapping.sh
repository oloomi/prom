#!/bin/bash
# Running SAMTools for BAM file generation

sam_file="corrected-mappings-mtb-10-runs-ref-base-700-1000"

# Creating BAM files for corrected mtb read mapping
samtools view -bS $sam_file.sam -o $sam_file.bam
samtools sort $sam_file.bam $sam_file.sorted
samtools index $sam_file.sorted.bam

echo "\n=== Done ===\n"
