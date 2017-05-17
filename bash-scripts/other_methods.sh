#!/bin/bash
# Running other multi-mapping resolution methods

reference="../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna"
alignments="../read-mapping/mtb-whole-genome-mutated/mtb-wg-mutated-se-mapping-report-all"
outfile="../read-mapping/mtb-whole-genome-mutated/corrected-other-3mis"

#reference="../data/genomes/mtb-genome-extract.fna"
#alignments="../read-mapping/mtb-mutated-long-repeats/mtb-mutated-se-mapping-report-all"
#outfile="../read-mapping/mtb-mutated-long-repeats/corrected-other-3mis"

export PATH=$PATH:/mnt/e/Codes/other-methods/mmr/mmr

samtools sort -n -o $alignments-id-sorted.bam $alignments.bam
#mmr -o $outfile.mmr.bam -f -b -R 150 $alignments.id-sorted.bam
#mmr -o $outfile.mmr.bam -F 0 -b -R 150 $alignments.id-sorted.bam

STARTTIME=$(date +%s)

mmr -o $outfile-mmr.bam -F 3 -b -R 150 $alignments-id-sorted.bam

ENDTIME=$(date +%s)
echo "MMR running time $(($ENDTIME - $STARTTIME)) seconds\n"

sam_file=$outfile-mmr
samtools sort $sam_file.bam -o $sam_file-sorted.bam
samtools index $sam_file-sorted.bam

echo "\n=== Trying Other Multimapping Resolution Tools Done ===\n"
