#!/bin/bash
# Running SAMTools for BAM file generation

#sam_file="../read-mapping/toy-genome-mutated/corrected-toy-wg-mutated-se-mapping-simple-bayesian"
sam_file="../read-mapping/toy-genome-mutated/corrected-toy-wg-mutated-se-mapping-remu"
#sam_file="../read-mapping/toy-genome-mutated-middle/corrected-toy-wg-mutated-middle-se-mapping-prom"
#sam_file="../read-mapping/toy-genome-mutated-middle/corrected-toy-wg-mutated-middle-se-mapping-remu"
#sam_file="../read-mapping/mtb-whole-genome-mutated-70-140/corrected-mtb-wg-mutated-se-mapping-remu"
#sam_file="../read-mapping/mtb-whole-genome-mutated-100-140/corrected-mtb-wg-mutated-se-mapping-remu-25-pmu"
#sam_file="../read-mapping/mtb-whole-genome-mutated-100-140/simple-bayesian-mtb-wg-mutated-se-mapping-25"
#sam_file="../read-mapping/mtb-whole-genome-mutated-70-140/mtb-wg-mutated-se-mapping-report-all-unique"
#sam_file="../read-mapping/mtb-whole-genome-mutated/corrected-mtb-wg-mutated-se-mapping"
#sam_file="../read-mapping/ot-whole-genome-mutated-70-140/corrected-ot-wg-mutated-se-mapping-remu"
#sam_file="../read-mapping/ot-whole-genome-mutated-70-140/corrected-ot-wg-mutated-se-mapping-remu-25-pmu"

# Creating BAM files for corrected mtb read mapping
samtools view -bS $sam_file.sam -o $sam_file.bam
samtools sort $sam_file.bam -o $sam_file-sorted.bam
samtools index $sam_file-sorted.bam

echo "\n=== Done ===\n"
