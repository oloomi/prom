# bayesian_update("/home/mohammad/pneumoniae/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated.fna",
#                 "/home/mohammad/pneumoniae/read-mapping/kt-kpnih1-back-mutated/kt-kpnih1-bm-report-all.sam",
#                 "/home/mohammad/pneumoniae/read-mapping/kt-kpnih1-back-mutated/remu-kt-kpnih1-bm-report-all.sam")

# bayesian_update("./data/genomes/toy-genome.fna",
#                 "./read-mapping/toy-genome-mutated/toy-wg-mutated-se-mapping-report-all.sam",
#                 "./read-mapping/toy-genome-mutated/corrected-toy-wg-mutated-se-mapping-remu.sam")

# bayesian_update("./data/genomes/toy-genome.fna",
#                 "./read-mapping/toy-genome-mutated-middle/toy-wg-mutated-middle-se-mapping-report-all.sam",
#                 "./read-mapping/toy-genome-mutated-middle/corrected-toy-wg-mutated-middle-se-mapping-remu.sam")

# bayesian_update("./data/genomes/Orientia_tsutsugamushi_Ikeda_uid58869/NC_010793.fna",
#                 "./read-mapping/ot-whole-genome-mutated-70-140/ot-wg-mutated-se-mapping-report-all.sam",
#                 "./read-mapping/ot-whole-genome-mutated-70-140/corrected-ot-wg-mutated-se-mapping-remu-25-pmu.sam")

# bayesian_update("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna",
#                 "./read-mapping/mtb-whole-genome-mutated-70-140/mtb-wg-mutated-se-mapping-report-all.sam",
#                 "./read-mapping/mtb-whole-genome-mutated-70-140/corrected-mtb-wg-mutated-se-mapping-remu.sam")

# bayesian_update("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna",
#                 "./read-mapping/mtb-whole-genome-mutated-100-140/mtb-wg-mutated-se-mapping-report-all.sam",
#                 "./read-mapping/mtb-whole-genome-mutated-100-140/corrected-mtb-wg-mutated-se-mapping-remu-25-pmu.sam")

# bayesian_update("./data/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full.fna",
#                     "./read-mapping/kp-kpnih1-back-mutated-full-real/kp-back-mutated-full-mapping-report-all.sam",
#                     "./read-mapping/kp-kpnih1-back-mutated-full-real/kp-back-mutated-full-mapping-remu.sam")

# bayesian_update("./data/genomes/MTB-H37Rv-back-mutated-full.fna",
#                 "./read-mapping/mtb-h37rv-back-mutated/mtb-h37rv-back-mutated-mapping-report-all.sam",
#                 "./read-mapping/mtb-h37rv-back-mutated/mtb-h37rv-back-mutated-mapping-remu.sam")


# file_path = "./read-mapping/ot-whole-genome-mutated-70-140/"
#
# vcf_files_names = [["Bowtie2 best-match", "ot-wg-mutated-se-mapping-best-match-sorted"],
#                    ["Bowtie2 report-all", "ot-wg-mutated-se-mapping-report-all-sorted"],
#                    ["MMR", "corrected-other-3mis-mmr-sorted"],
#                    ["REMU", "corrected-ot-wg-mutated-se-mapping-remu-sorted"],
#                    ["REMU-pmu", "corrected-ot-wg-mutated-se-mapping-remu-25-pmu-sorted"]]
#
# evaluation_results = open("./results/variants-comparison-OT-wg-70-140-pmu25.txt", 'w')

# comparison_output = \
    #     compare_variants("/mnt/e/Codes/bayesian-update/data/genomes/ot-whole-genome-mutated-70-140-mutations.txt",
    #                      vcf_files)

# find_unique_reads("./read-mapping/mtb-whole-genome-mutated-70-140/mtb-wg-mutated-se-mapping-report-all.sam")

# remu.py

    # process = psutil.Process(os.getpid())
    # mem = process.memory_info()[0] / float(2 ** 20)
    # print("Memory usage:", round(mem, 2), "MB")


# -r "./data/genomes/toy-genome.fna" -i "./read-mapping/toy-genome-mutated/toy-wg-mutated-se-mapping-report-all.sam" -o "./read-mapping/toy-genome-mutated/test.sam"

# python3 remu.py -r "../data/NC_000962.fna" -i "../data/toy-report-all.sam" -o "../data/remu-toy-report-all.sam"

#reads_genome="../data/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full.fna"
#reads_genome="../data/genomes/toy-genome-mutated-middle.fna"
#reads_genome="../data/genomes/toy-genome-mutated-middle.fna"
#reads_genome="../data/genomes/ot-whole-genome-mutated-70-140.fna"
#reads_genome="../data/genomes/mtb-whole-genome-mutated-70-140.fna"
#reads_genome="../data/genomes/mtb-whole-genome-mutated-100-140-half.fna"
#reads_genome="../data/genomes/mtb-whole-genome-mutated.fna"
#reads_genome="../data/genomes/mtb-genome-extract.fna"
#reads_genome="../data/genomes/mtb-genome-extract-mutated-long-repeats.fna"

#out_dir="../read-mapping/toy-genome-mutated/"
#file_prefix="toy-wg-mutated-se"
#out_dir="../read-mapping/toy-genome-mutated-middle/"
#file_prefix="toy-wg-mutated-middle-se"
#out_dir="../read-mapping/ot-whole-genome-mutated-70-140/"
#file_prefix="ot-wg-mutated-se"
#out_dir="../read-mapping/mtb-whole-genome-mutated-100-140/"
#file_prefix="mtb-wg-mutated-se"
#out_dir="../read-mapping/mtb-whole-genome-mutated-70-140/"
#out_dir="../read-mapping/mtb-whole-genome-mutated/"
#file_prefix="mtb-wg-mutated-se"
#out_dir="../read-mapping/mtb-mutated-long-repeats/"
#file_prefix="mtb-mutated-se"
#out_dir="../read-mapping/mtb-normal/"
#file_prefix="mtb-normal-se"

#reference="../data/genomes/Orientia_tsutsugamushi_Ikeda_uid58869/NC_010793.fna"
#alignments="../read-mapping/ot-whole-genome-mutated-70-140/ot-wg-mutated-se-mapping-report-all"
#outfile="../read-mapping/ot-whole-genome-mutated-70-140/corrected-other-3mis"

#reference="../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna"
#alignments="../read-mapping/mtb-whole-genome-mutated-100-140/mtb-wg-mutated-se-mapping-report-all"
#outfile="../read-mapping/mtb-whole-genome-mutated-100-140/corrected-other-3mis"
#alignments="../read-mapping/mtb-whole-genome-mutated-70-140/mtb-wg-mutated-se-mapping-report-all"
#outfile="../read-mapping/mtb-whole-genome-mutated-70-140/corrected-other-3mis"

#reference="../data/genomes/mtb-genome-extract.fna"
#alignments="../read-mapping/mtb-mutated-long-repeats/mtb-mutated-se-mapping-report-all"
#outfile="../read-mapping/mtb-mutated-long-repeats/corrected-other-3mis"

#alignment_files="mtb-h37rv-back-mutated-mapping-mmr-sorted"

#reference="../data/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full.fna"
#file_path="../read-mapping/kp-kpnih1-back-mutated-full-real/"

#alignment_files="kp-back-mutated-full-mapping-best-match-sorted
#kp-back-mutated-full-mapping-report-all-sorted
#kp-back-mutated-full-mapping-prom-sorted
#kp-back-mutated-full-mapping-remu-sorted"

#alignment_files="kp-back-mutated-full-mapping-prom-sorted"

#reference="../data/genomes/Orientia_tsutsugamushi_Ikeda_uid58869/NC_010793.fna"
#file_path="../read-mapping/ot-whole-genome-mutated-70-140/"

#alignment_files="corrected-ot-wg-mutated-se-mapping-remu-sorted"

#alignment_files="corrected-other-3mis-mmr-sorted"

#alignment_files="ot-wg-mutated-se-sorted
#ot-wg-mutated-se-mapping-best-match-sorted
#ot-wg-mutated-se-mapping-report-all-sorted
#corrected-other-3mis-mmr-sorted
#corrected-ot-wg-mutated-se-mapping-remu-sorted"

#reference="../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna"
#file_path="../read-mapping/mtb-whole-genome-mutated-70-140/"
#file_path="../read-mapping/mtb-whole-genome-mutated-100-140/"

#alignment_files="corrected-ot-wg-mutated-se-mapping-remu-25-pmu-sorted"

#alignment_files="simple-bayesian-mtb-wg-mutated-se-mapping-25-sorted"

#alignment_files="mtb-wg-mutated-se-sorted
#mtb-wg-mutated-se-mapping-best-match-sorted
#mtb-wg-mutated-se-mapping-report-all-sorted
#corrected-other-3mis-mmr-sorted
#corrected-mtb-wg-mutated-se-mapping-remu-25-pmu-sorted"


#alignments="../read-mapping/mtb-mutated-long-repeats/corrected-mappings-mtb-mutated-700-100-1-10runs.sorted"

#bcftools mpileup -f $reference $alignments.bam | bcftools call -mv --ploidy 1 -P 1.1e-1 -o $alignments-variants.vcf

#reference="../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna"
#file_path="../read-mapping/mtb-h37rv-back-mutated/"
#alignment_files="mtb-h37rv-original-mapping-best-match-sorted"

#	bcftools mpileup -f $reference $file_path$file.bam | bcftools call -mv --ploidy 1 -P 1.1e-1 -o $file_path$file-variants-mv.vcf
#	bcftools mpileup -f $reference $file_path$file.bam | bcftools call -cv --ploidy 1 -p 0.5 -o $file_path$file-variants-consensus-p0.5.vcf
#	freebayes -f $reference -p 1 $file_path$file.bam >$file_path$file-variants-freebayes.vcf

# variant_caller_lst = [("Freebayes", "freebayes"), ("FreebayesMin", "freebayes-min"),
    #                       ("BCFtools p 0.5", "consensus-p0.5"), ("BCFtools mv", "mv")]

# vcf_files_names = [["Bowtie2 best-match", "mtb-wg-mutated-se-mapping-best-match-sorted"],
#                        ["Bowtie2 report-all", "mtb-wg-mutated-se-mapping-report-all-sorted"],
#                        ["MMR", "corrected-other-3mis-mmr-sorted"],
#                        ["PROM", "simple-bayesian-mtb-wg-mutated-se-mapping-25-sorted"],
#                        ["REMU", "corrected-mtb-wg-mutated-se-mapping-remu-25-pmu-sorted"]]