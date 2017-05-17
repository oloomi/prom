#!/bin/bash
# Running BCFtools for variant calling

#reference="../data/genomes/mtb-genome-extract.fna"
reference="../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna"

#file_path="../read-mapping/mtb-mutated-long-repeats/"
file_path="../read-mapping/mtb-whole-genome-mutated/"

#alignment_files="mtb-wg-mutated-se-sorted
#mtb-wg-mutated-se-mapping-best-match-sorted
#mtb--wg-mutated-se-mapping-report-all-sorted
#corrected-other.mmr.sorted
#corrected-other-3mis.mmr.sorted
#corrected-other-best.mmr.sorted
#corrected-mtb-wg-mutated-se-mapping.sorted"

alignment_files="mtb-wg-mutated-se-sorted
mtb-wg-mutated-se-mapping-best-match-sorted
mtb-wg-mutated-se-mapping-report-all-sorted
corrected-other-3mis-mmr-sorted
corrected-mtb-wg-mutated-se-mapping-sorted"


#alignments="../read-mapping/mtb-mutated-long-repeats/corrected-mappings-mtb-mutated-700-100-1-10runs.sorted"

#bcftools mpileup -f $reference $alignments.bam | bcftools call -mv --ploidy 1 -P 1.1e-1 -o $alignments-variants.vcf

for file in $alignment_files
do
	echo "\nVariant calling for: $file\n"
#	bcftools mpileup -f $reference $file_path$file.bam | bcftools call -mv --ploidy 1 -P 1.1e-1 -o $file_path$file-variants.vcf
#	bcftools mpileup -f $reference $file_path$file.bam | bcftools call -cv --ploidy 1 -p 0.99 -o $file_path$file-variants-consensus-p0.99.vcf
	bcftools mpileup -f $reference $file_path$file.bam | bcftools call -cv --ploidy 1 -p 0.5 -o $file_path$file-variants-consensus-p0.5.vcf
#	freebayes -f $reference -p 1 $file_path$file.bam >$file_path$file-variants-freebayes.vcf
done

echo "\n=== Done ===\n"
