#!/bin/bash
# Running SAMTools for BAM file generation

mtb_genome="./mtb-genome-extract.fna"

# Creating BAM files for MTB ArtIllumina benchmark
samtools view -bS mtb-single-end.sam -o mtb-single-end.bam
samtools sort mtb-single-end.bam mtb-single-end-sorted
samtools index mtb-single-end-sorted.bam

# Creating BAM files for MTB best-match
samtools view -bS mtb-single-end-mapping-best-match.sam -o mtb-single-end-mapping-best-match.bam
samtools sort mtb-single-end-mapping-best-match.bam mtb-single-end-mapping-best-match-sorted
samtools index mtb-single-end-mapping-best-match-sorted.bam

# Creating BAM files for MTB report-all
samtools view -bS mtb-single-end-mapping-report-all.sam -o mtb-single-end-mapping-report-all.bam
samtools sort mtb-single-end-mapping-report-all.bam mtb-single-end-mapping-report-all-sorted
samtools index mtb-single-end-mapping-report-all-sorted.bam


echo "\n=== Done ===\n"
