#!/usr/bin/python3

import sys
import copy
from genome_util import *
from vcf_file import *


def create_mutated_gnome():
    """
    Run this function in the experiment main folder eg. supermax-100-140
    """
    mutate_genome_repeats("../../genome-ref/ref-genome.fna",
                          "../../genome-ref/repeats/supermax-repeats.txt",
                          "./genome-mutated/mutated-genome.fna",
                          "./genome-mutated/mutations.txt")
    return True


def variant_evaluation():
    variant_caller_lst = [("Freebayes", "freebayes")]

    file_path = "./variants/"
    # file_path = "./read-mapping/mtb-whole-genome-mutated-100-140/"

    vcf_files_names = [["Bowtie2 best-match", "bowtie-mapping-best-match-sorted"],
                       ["Bowtie2 report-all", "bowtie-mapping-report-all-sorted"],
                       ["Bowtie2 + MMR", "bowtie-mmr-sorted"],
                       ["Bowtie2 + REMU", "bowtie-remu-sorted"],
                       ["BWA best-match", "bwa-mapping-best-match-sorted"],
                       ["BWA report-all", "bwa-mapping-report-all-sorted"],
                       ["BWA + MMR", "bwa-mmr-sorted"],
                       ["BWA + REMU", "bwa-remu-sorted"]]

    evaluation_results = open("./results/variants-comparison-freebayes.txt", 'w')

    for variant_caller in variant_caller_lst:
        vcf_files = copy.deepcopy(vcf_files_names)
        for i in range(len(vcf_files)):
            vcf_files[i][1] = file_path + vcf_files[i][1] + "-variants-{}.vcf".format(variant_caller[1])

        comparison_output = \
            compare_variants(
                "./genome-mutated/mutations.txt",
                vcf_files)

        evaluation_results.write("{}\n".format(variant_caller[0]))
        evaluation_results.write(comparison_output)

    evaluation_results.close()


if __name__ == "__main__":
    if sys.argv[1] == '-c':
        create_mutated_gnome()
    elif sys.argv[1] == '-e':
        variant_evaluation()
    else:
        print("Invalid argument!")