#!/usr/bin/env bash

# Run the following command in the experiment folder after all stages are finished
# sh bwa-mmr-fix.sh [-s -r] [read_len]

out_dir="./mappings/bwa/"
file_prefix="bwa"

# Completing sequence and quality score fields for secondary alignments
echo "Completing BWA fields for secondary alignments...\n"
script_path="$( cd "$(dirname "$0")" ; pwd -P )/"
python3 ${script_path}"bwa-complete-secondary.py" $2 ${out_dir}${file_prefix}-mapping-report-all.sam ${out_dir}${file_prefix}-mapping-report-all-complete.sam

# Creating BAM files for report-all
samtools view -bS ${out_dir}${file_prefix}-mapping-report-all-complete.sam -o ${out_dir}${file_prefix}-mapping-report-all-complete.bam
samtools sort ${out_dir}${file_prefix}-mapping-report-all-complete.bam -o ${out_dir}${file_prefix}-mapping-report-all-complete-sorted.bam
samtools index ${out_dir}${file_prefix}-mapping-report-all-complete-sorted.bam

# MMR method
echo "Running MMR...\n"
alignments="./mappings/bwa/bwa-mapping-report-all-complete-sorted"
outfile="./mappings/bwa/bwa"

samtools sort -n ${alignments}.bam -o ${alignments}-id-sorted.bam

/usr/bin/time -v -o ${outfile}-mmr-time-log-2.txt mmr -o ${outfile}-mmr-fixed.bam -F 3 -b -R $2 ${alignments}-id-sorted.bam | tee ${outfile}-mmr-log-2.txt

sam_file=${outfile}-mmr-fixed
samtools sort ${sam_file}.bam -o ${sam_file}-sorted.bam
samtools index ${sam_file}-sorted.bam

# Variant calling
if [ "$1" = "-s" ]; then
  ref_genome="../../genome-ref/ref-genome.fna"
elif [ "$1" = "-r" ]; then
  ref_genome="./genome-mutated/mutated-genome.fna"
else
  echo "Invalid argument in variant-calling.sh!"
fi

file_path="./mappings/bwa/"
alignment_files="bwa-mmr-fixed-sorted"
out_path="./variants/"


for file in ${alignment_files}
do
	echo "\nVariant calling for: $file\n"
#	-m 0 since BWA gives mapping quality zero to multireads
	freebayes -f ${ref_genome} -p 1 -F 0.9 -m 0 ${file_path}${file}.bam >${out_path}bwa-mmr-sorted-variants-freebayes.vcf
done

echo "\n=== BWA+MMR fix completed! ===\n"

# Re-running evaluation
python3 ${script_path}"evaluation.py" -e