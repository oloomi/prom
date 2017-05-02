#!/bin/bash
# Running BCFtools for variant calling

reference="../data/genomes/mtb-genome-extract.fna"

file_path="../read-mapping/mtb-mutated-long-repeats/"

alignment_files="mtb-mutated-se-sorted
mtb-mutated-se-mapping-best-match-sorted
mtb-mutated-se-mapping-report-all-sorted
corrected-other.mmr.sorted
corrected-other-3mis.mmr.sorted
corrected-other-best.mmr.sorted
corrected-mappings-mtb-mutated-700-100-1-10runs.sorted"

#alignments="../read-mapping/mtb-mutated-long-repeats/corrected-mappings-mtb-mutated-700-100-1-10runs.sorted"

#bcftools mpileup -f $reference $alignments.bam | bcftools call -mv --ploidy 1 -P 1.1e-1 -o $alignments-variants.vcf

for file in $alignment_files
do
#	bcftools mpileup -f $reference $file_path$file.bam | bcftools call -mv --ploidy 1 -P 1.1e-1 -o $file_path$file-variants.vcf
#	bcftools mpileup -f $reference $file_path$file.bam | bcftools call -cv --ploidy 1 -p 0.99 -o $file_path$file-variants-consensus-p0.99.vcf
	freebayes -f $reference -p 1 $file_path$file.bam >$file_path$file-variants-freebayes.vcf 
done

echo "\n=== Done ===\n"
