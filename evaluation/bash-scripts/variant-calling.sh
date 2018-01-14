#!/bin/bash

reference="../../genome-ref/ref-genome.fna"
out_path="./variants/"

file_path="./mappings/bowtie/"
alignment_files="bowtie-mapping-best-match-sorted
bowtie-mapping-report-all-sorted
bowtie-mmr-sorted
bowtie-remu-sorted"


for file in $alignment_files
do
	echo "\nVariant calling for: $file\n"
	freebayes -f $reference -p 1 -F 0.9 $file_path$file.bam >$out_path$file-variants-freebayes.vcf
done

file_path="./mappings/bwa/"
alignment_files="bwa-mapping-best-match-sorted
bwa-mapping-report-all-sorted
bwa-mmr-sorted
bwa-remu-sorted"


for file in $alignment_files
do
	echo "\nVariant calling for: $file\n"
#	-m 0 since BWA gives mapping quality zero to multireads
	freebayes -f $reference -p 1 -F 0.9 -m 0 $file_path$file.bam >$out_path$file-variants-freebayes.vcf
done

echo "\n=== Variant Calling Done ===\n"
