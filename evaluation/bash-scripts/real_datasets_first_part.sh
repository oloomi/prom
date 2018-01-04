# the reference genome sequence used for read mapping

original_genome="../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna"
#ref_genome="../data/genomes/MTB-H37Rv-back-mutated-full.fna"

reads_1="../read-datasets/MTB-H37Rv/fastq/SRR2818101_pass_1_trim_250.fastq"

out_dir="../read-mapping/mtb-h37rv-back-mutated/"
#file_prefix="mtb-h37rv-back-mutated"
file_prefix="mtb-h37rv-original"

#ref_genome="../data/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full.fna"
#ref_genome="../data/genomes/Klebsiella_pneumoniae_KPNIH1/Klebsiella_pneumoniae_KPNIH1.fna"

#reads_1="../read-datasets/KPNIH1/SRR1505904_pass_1.fastq"
#reads_2="../read-datasets/KPNIH1/SRR1505904_pass_1.fastq"

#out_dir="../read-mapping/kp-kpnih1-back-mutated-full-real/"
#file_prefix="kp-back-mutated-full"

# ==== Read mapping ====

# Build index for original genome
bowtie2-build $original_genome $out_dir"original-genome-index"

# Single mapping, best-match for MTB original genome
bowtie2 -x $out_dir"original-genome-index" -U $reads_1,$reads_2 -S $out_dir$file_prefix-mapping-best-match.sam 2> $out_dir$file_prefix-mapping-best-match-log.txt

# Build index for reference genome
#bowtie2-build $ref_genome $out_dir"genome-index"

# Single mapping, best-match for MTB
#bowtie2 -x $out_dir"genome-index" -U $reads_1,$reads_2 -S $out_dir$file_prefix-mapping-best-match.sam 2> $out_dir$file_prefix-mapping-best-match-log.txt

# Single mapping, report-all for MTB
#bowtie2 -a -x $out_dir"genome-index" -U $reads_1,$reads_2 -S $out_dir$file_prefix-mapping-report-all.sam 2> $out_dir$file_prefix-mapping-report-all-log.txt

echo "\n=== Read mapping DONE! ===\n"

samtools view -bS $out_dir$file_prefix-mapping-best-match.sam | samtools sort - -o $out_dir$file_prefix-mapping-best-match-sorted.bam
samtools index $out_dir$file_prefix-mapping-best-match-sorted.bam


