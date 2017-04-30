#!/bin/bash
# Running BCFtools for variant calling

reference="../data/genomes/mtb-genome-extract.fna"
alignments="../read-mapping/mtb-mutated/corrected-mappings-mtb-mutated-700-100-1-10runs-fs.sorted"

bcftools mpileup -f $reference $alignments.bam | bcftools call -mv --ploidy 1 -P 1.1e-1 -o $alignments-variants.vcf

echo "\n=== Done ===\n"
