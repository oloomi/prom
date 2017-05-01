#!/bin/bash
# Running SAMTools for BAM file generation

#sam_file="../read-mapping/mtb-mutated-long-repeats/corrected-mappings-mtb-mutated-700-100-1-10runs-max"
#sam_file="../read-mapping/mtb-mutated-long-repeats/mtb-mutated-se-mapping-report-all-unique"
sam_file="../read-mapping/mtb-mutated-long-repeats/corrected-other-3mis.mmr"

# Creating BAM files for corrected mtb read mapping
#samtools view -bS $sam_file.sam -o $sam_file.bam
samtools sort $sam_file.bam -o $sam_file.sorted.bam
samtools index $sam_file.sorted.bam

echo "\n=== Done ===\n"
