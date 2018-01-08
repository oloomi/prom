import copy
from evaluation.vcf_file import *


variant_caller_lst = [("Freebayes", "freebayes"), ("FreebayesMin", "freebayes-min"),
                      ("BCFtools p 0.5", "consensus-p0.5"), ("BCFtools mv", "mv")]

file_path = "./read-mapping/mtb-whole-genome-mutated-100-140/"

vcf_files_names = [["Bowtie2 best-match", "mtb-wg-mutated-se-mapping-best-match-sorted"],
                   ["Bowtie2 report-all", "mtb-wg-mutated-se-mapping-report-all-sorted"],
                   ["MMR", "corrected-other-3mis-mmr-sorted"],
                   ["PROM", "simple-bayesian-mtb-wg-mutated-se-mapping-25-sorted"],
                   ["REMU","corrected-mtb-wg-mutated-se-mapping-remu-25-pmu-sorted"]]

evaluation_results = open("./results/variants-comparison-MTB-wg-100-140-freebayes-min.txt", 'w')

for variant_caller in variant_caller_lst:
    vcf_files = copy.deepcopy(vcf_files_names)
    for i in range(len(vcf_files)):
        vcf_files[i][1] = file_path + vcf_files[i][1] + "-variants-{}.vcf".format(variant_caller[1])

    comparison_output = \
        compare_variants("/mnt/e/Codes/bayesian-update/data/genomes/mtb-whole-genome-mutated-100-140-half-mutations.txt",
                         vcf_files)

    evaluation_results.write("{}\n".format(variant_caller[0]))
    evaluation_results.write(comparison_output)

evaluation_results.close()
