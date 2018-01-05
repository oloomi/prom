import copy
import timeit

from evaluation.vcf_file import *


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


variant_caller_lst = [("Freebayes", "freebayes"), ("FreebayesMin", "freebayes-min"), ("BCFtools p 0.5", "consensus-p0.5"), ("BCFtools mv", "mv")]

file_path = "./read-mapping/mtb-whole-genome-mutated-100-140/"

vcf_files_names = [["Bowtie2 best-match", "mtb-wg-mutated-se-mapping-best-match-sorted"],
                   ["Bowtie2 report-all", "mtb-wg-mutated-se-mapping-report-all-sorted"],
                   ["MMR", "corrected-other-3mis-mmr-sorted"],
                   ["PROM", "simple-bayesian-mtb-wg-mutated-se-mapping-25-sorted"],
                   ["REMU","corrected-mtb-wg-mutated-se-mapping-remu-25-pmu-sorted"]]

evaluation_results = open("./results/variants-comparison-MTB-wg-100-140-freebayes-min.txt", 'w')

# file_path = "./read-mapping/ot-whole-genome-mutated-70-140/"
#
# vcf_files_names = [["Bowtie2 best-match", "ot-wg-mutated-se-mapping-best-match-sorted"],
#                    ["Bowtie2 report-all", "ot-wg-mutated-se-mapping-report-all-sorted"],
#                    ["MMR", "corrected-other-3mis-mmr-sorted"],
#                    ["REMU", "corrected-ot-wg-mutated-se-mapping-remu-sorted"],
#                    ["REMU-pmu", "corrected-ot-wg-mutated-se-mapping-remu-25-pmu-sorted"]]
#
# evaluation_results = open("./results/variants-comparison-OT-wg-70-140-pmu25.txt", 'w')

for variant_caller in variant_caller_lst:
    vcf_files = copy.deepcopy(vcf_files_names)
    for i in range(len(vcf_files)):
        vcf_files[i][1] = file_path + vcf_files[i][1] + "-variants-{}.vcf".format(variant_caller[1])

    comparison_output = \
        compare_variants("/mnt/e/Codes/bayesian-update/data/genomes/mtb-whole-genome-mutated-100-140-half-mutations.txt",
                         vcf_files)
    # comparison_output = \
    #     compare_variants("/mnt/e/Codes/bayesian-update/data/genomes/ot-whole-genome-mutated-70-140-mutations.txt",
    #                      vcf_files)
    evaluation_results.write("{}\n".format(variant_caller[0]))
    evaluation_results.write(comparison_output)

evaluation_results.close()

# find_unique_reads("./read-mapping/mtb-whole-genome-mutated-70-140/mtb-wg-mutated-se-mapping-report-all.sam")