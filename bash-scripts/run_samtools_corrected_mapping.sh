#!/bin/bash
# Running SAMTools for BAM file generation

#sam_file="../read-mapping/ot-whole-genome-mutated-70-140/corrected-ot-wg-mutated-se-mapping-filter"
sam_file="../read-mapping/mtb-whole-genome-mutated-70-140/corrected-mtb-wg-mutated-se-mapping-filter"
#sam_file="../read-mapping/mtb-whole-genome-mutated-70-140/mtb-wg-mutated-se-mapping-report-all-unique"
#sam_file="../read-mapping/mtb-whole-genome-mutated/corrected-mtb-wg-mutated-se-mapping"

# Creating BAM files for corrected mtb read mapping
samtools view -bS $sam_file.sam -o $sam_file.bam
samtools sort $sam_file.bam -o $sam_file-sorted.bam
samtools index $sam_file-sorted.bam

echo "\n=== Done ===\n"
