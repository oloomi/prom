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