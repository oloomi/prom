from sam_file import *
from select_mapping import *
from calc_likelihood import *
from vcf_file import *
import random
import copy
import timeit


def bayesian_resolution(ref_genome_file, sam_file, output_file):
    """
    Assigning a multi-read to a mapping location using Bayesian updating
    :param ref_genome_file:
    :param sam_file:
    :param output_file:
    :return:
    """

    # 0. Reading reference genome FASTA file and mapping SAM file
    genome_header, genome_seq = read_genome(ref_genome_file)
    reads_dict, read_len = read_sam_file(sam_file, genome_seq)

    # Estimated average depth of coverage
    coverage = int(len(reads_dict) * read_len / len(genome_seq))
    print("Estimated average depth of coverage according to mapped reads: {}".format(coverage))

    # 1. Finding initial counts
    initial_base_counts = initial_counts(genome_seq, coverage)

    multi_reads_final_location = defaultdict(int)

    # Updating the prior (initial counts) for uniquely mapped reads
    # New: and also filtering may make a multi-read a unique read
    unique_reads = defaultdict(int)
    initially_resolved_multireads = defaultdict(int)
    random.seed(12)

    for read_id, mappings in reads_dict.items():
        if len(mappings) == 1:
            unique_reads[read_id] = mappings[0][0]
            mapping_start_pos = mappings[0][0] - 1
            read_seq = mappings[0][3]

            update_counts(initial_base_counts, [0.99, mapping_start_pos, read_seq], coverage)     # 0.99 won't be used

            # for pos_in_read, base in enumerate(read_seq):
            #     # We find the base counts for that position in reference genome and
            #     # we update it by incrementing the count for the base in the read
            #     initial_base_counts[mapping_start_pos + pos_in_read][base_index[base]] += 1
        else:
            mappings_filtered = filter_alignments(mappings, 3)
            # By removing rubbish alignments, it has turned to a unique mapping
            if len(mappings_filtered) == 1:
                multi_reads_final_location[read_id] = mappings_filtered[0][0]
                initially_resolved_multireads[read_id] = mappings_filtered[0][0]
                # We treat it like a unique read and update counts
                mapping_start_pos = mappings_filtered[0][0] - 1
                read_seq = mappings_filtered[0][3]

                update_counts(initial_base_counts, [0.99, mapping_start_pos, read_seq], coverage)

                # for pos_in_read, base in enumerate(read_seq):
                #     initial_base_counts[mapping_start_pos + pos_in_read][base_index[base]] += 1
            else:
                # It is a multi-mapping that needs to be resolved by Bayesian updating
                # However, rubbish mapping locations should be removed
                reads_dict[read_id] = mappings_filtered

    # Writing unique reads in SAM file
    # unique_reads_write_sam(sam_file, unique_reads, initially_resolved_multireads)
    # print("Writing unique reads done!")
    # return True

    # Removing uniquely mapped reads
    for read_id in unique_reads:
        del reads_dict[read_id]

    print("Number of multireads: {}".format(len(reads_dict)))

    # Removing initially resolved multireads
    for read_id in initially_resolved_multireads:
        del reads_dict[read_id]

    print("Number of initially resolved multireads: {}".format(len(initially_resolved_multireads)))

    with open("multireads.txt", "w") as multireads_file:
        for key, value in reads_dict.items():
            mdz_s = [lst[2] for lst in value]
            multireads_file.write("{}\t{}\n".format(key, mdz_s))

    multi_reads = sorted(reads_dict.keys())
    print("Number of multi-reads requiring resolution:", len(multi_reads))

    for run_number in range(1):
        print("Run # {} ...".format(run_number))
        random.seed(run_number)
        base_counts = copy.deepcopy(initial_base_counts)

        for read_id in reads_dict:

            # [(probability, position, read_seq), ...] where read_seq is needed for updating counts later on
            mapping_probs = []

            # For each of its mapping location, we calculate the posterior probability
            for mapping in reads_dict[read_id]:
                (mapping_start_pos, read_seq) = (mapping[0] - 1, mapping[3])
                # Find the likelihood
                mapping_probs.append([calc_log_mapping_prob(base_counts, mapping_start_pos, read_seq),
                                      mapping_start_pos, read_seq])

            # Selecting one location statistically
            # After select_mapping, the list `mapping_probs` is modified and contains probabilities instead of log-probs
            selected_mapping = best_mapping(mapping_probs)

            multi_reads_final_location[read_id] = selected_mapping[1] + 1

            # Updating base counts for selected location
            # update_counts(base_counts, selected_mapping)

    # Writing final results to a SAM file
    write_sam_file(multi_reads_final_location, sam_file, output_file)

    return True


# bayesian_resolution("./data/genomes/mtb-genome-extract.fna",
#                 "./read-mapping/mtb-normal/mtb-normal-se-mapping-report-all.sam",
#                 "./read-mapping/mtb-normal/corrected-mappings-mtb-normal-700-100-5.sam")

# bayesian_resolution("./data/genomes/mtb-genome-extract-mutated.fna",
#                 "./read-mapping/mtb-mutated/mtb-mutated-se-mapping-report-all.sam",
#                 "./read-mapping/mtb-mutated/corrected-mappings-mtb-mutated-700-100-1-10runs-fs.sam")

# bayesian_resolution("./data/genomes/mtb-genome-extract.fna",
#                 "./read-mapping/mtb-mutated-long-repeats/mtb-mutated-se-mapping-report-all.sam",
#                 "./read-mapping/mtb-mutated-long-repeats/corrected-mappings-mtb-mutated-700-100-1-10runs-max.sam")

phase = 2

if phase == 1:
    start_time = timeit.default_timer()

    # bayesian_resolution("./data/genomes/toy-genome.fna",
    #                 "./read-mapping/toy-genome-mutated/toy-wg-mutated-se-mapping-report-all.sam",
    #                 "./read-mapping/toy-genome-mutated/corrected-toy-wg-mutated-se-mapping-simple-bayesian.sam")

    # bayesian_resolution("./data/genomes/toy-genome.fna",
    #                 "./read-mapping/toy-genome-mutated-middle/toy-wg-mutated-middle-se-mapping-report-all.sam",
    #                 "./read-mapping/toy-genome-mutated-middle/corrected-toy-wg-mutated-middle-se-mapping-prom.sam")

    # bayesian_resolution("./data/genomes/Orientia_tsutsugamushi_Ikeda_uid58869/NC_010793.fna",
    #                 "./read-mapping/ot-whole-genome-mutated-70-140/ot-wg-mutated-se-mapping-report-all.sam",
    #                 "./read-mapping/ot-whole-genome-mutated-70-140/corrected-ot-wg-mutated-se-mapping-remu.sam")

    # bayesian_resolution("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna",
    #                 "./read-mapping/mtb-whole-genome-mutated-70-140/mtb-wg-mutated-se-mapping-report-all.sam",
    #                 "./read-mapping/mtb-whole-genome-mutated-70-140/corrected-mtb-wg-mutated-se-mapping-remu.sam")
    #
    # bayesian_resolution("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna",
    #                 "./read-mapping/mtb-whole-genome-mutated-100-140/mtb-wg-mutated-se-mapping-report-all.sam",
    #                 "./read-mapping/mtb-whole-genome-mutated-100-140/simple-bayesian-mtb-wg-mutated-se-mapping-25.sam")

    # bayesian_resolution("/home/mohammad/pneumoniae/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated.fna",
    #                 "/home/mohammad/pneumoniae/read-mapping/kt-kpnih1-back-mutated/kt-kpnih1-bm-report-all.sam",
    #                 "/home/mohammad/pneumoniae/read-mapping/kt-kpnih1-back-mutated/prom-kt-kpnih1-bm-report-all.sam")

    # bayesian_resolution("/home/mohammad/pneumoniae/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full.fna",
    #                     "/home/mohammad/pneumoniae/read-mapping/kt-kpnih1-back-mutated-full/kt-kpnih1-bm-report-all.sam",
    #                     "/home/mohammad/pneumoniae/read-mapping/kt-kpnih1-back-mutated-full/prom-kt-kpnih1-bm-report-all.sam")

    # bayesian_resolution("./data/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full.fna",
    #                     "./read-mapping/kp-kpnih1-back-mutated-full-real/kp-back-mutated-full-mapping-report-all.sam",
    #                     "./read-mapping/kp-kpnih1-back-mutated-full-real/kp-back-mutated-full-mapping-prom.sam")

    bayesian_resolution("./data/genomes/MTB-H37Rv-back-mutated-full.fna",
                        "./read-mapping/mtb-h37rv-back-mutated/mtb-h37rv-back-mutated-mapping-report-all.sam",
                        "./read-mapping/mtb-h37rv-back-mutated/mtb-h37rv-back-mutated-mapping-prom.sam")

    run_time = timeit.default_timer() - start_time

    print("Bayesian update running time: {} seconds = {} minutes".format(round(run_time, 2), round(run_time / 60, 2)))

# find_unique_reads("./read-mapping/mtb-mutated-long-repeats/mtb-mutated-se-mapping-report-all.sam")

elif phase == 2:

    variant_caller_lst = [("Freebayes", "freebayes"), ("FreebayesMin", "freebayes-min"), ("BCFtools p 0.5", "consensus-p0.5"), ("BCFtools mv", "mv")]

    file_path = "/mnt/e/Codes/bayesian-update/read-mapping/mtb-h37rv-back-mutated/"

    vcf_files_names = [["Bowtie2_best-match", "mtb-h37rv-back-mutated-mapping-best-match-sorted"],
                       ["Bowtie2_report-all", "mtb-h37rv-back-mutated-mapping-report-all-sorted"],
                       ["PROM", "mtb-h37rv-back-mutated-mapping-prom-sorted"],
                       ["REMU", "mtb-h37rv-back-mutated-mapping-remu-sorted"]
                       ]

    evaluation_results = open("./results/variants-comparison-mtb-h37rv-back-mutated-full-real-q20-new.txt", 'w')

    # variant_caller_lst = [("Freebayes", "freebayes"), ("BCFtools p 0.5", "consensus-p0.5"), ("BCFtools mv", "mv")]
    #
    # file_path = "/mnt/e/Codes/bayesian-update/read-mapping/kp-kpnih1-back-mutated-full-real/"
    #
    # vcf_files_names = [["Bowtie2_best-match", "kp-back-mutated-full-mapping-best-match-sorted"],
    #                    ["Bowtie2_report-all", "kp-back-mutated-full-mapping-report-all-sorted"],
    #                    ["PROM", "kp-back-mutated-full-mapping-prom-sorted"],
    #                    ["REMU", "kp-back-mutated-full-mapping-remu-sorted"]
    #                    ]
    #
    # evaluation_results = open("./results/variants-comparison-kp-kpnih1-back-mutated-full-real-new.txt", 'w')

    # variant_caller_lst = [("Freebayes", "freebayes"), ("BCFtools p 0.5", "consensus-p0.5"), ("BCFtools mv", "mv")]
    #
    # file_path = "./read-mapping/mtb-whole-genome-mutated-100-140/"
    #
    # vcf_files_names = [["Bowtie2 best-match", "mtb-wg-mutated-se-mapping-best-match-sorted"],
    #                    ["Bowtie2 report-all", "mtb-wg-mutated-se-mapping-report-all-sorted"],
    #                    ["MMR", "corrected-other-3mis-mmr-sorted"],
    #                    ["Bayesian", "simple-bayesian-mtb-wg-mutated-se-mapping-sorted"],
    #                    ["Bayesian-25", "simple-bayesian-mtb-wg-mutated-se-mapping-25-sorted"],
    #                    ["REMU", "corrected-mtb-wg-mutated-se-mapping-remu-sorted"],
    #                    ["REMU-25", "corrected-mtb-wg-mutated-se-mapping-remu-25-sorted"],
    #                    ["REMU-25-pmu", "corrected-mtb-wg-mutated-se-mapping-remu-25-pmu-sorted"]
    #                    ]
    #
    # evaluation_results = open("./results/variants-comparison-MTB-wg-100-140-simple-bayesian-25-pmu.txt", 'w')

    # file_path = "./read-mapping/ot-whole-genome-mutated-70-140/"
    #
    # vcf_files_names = [["Bowtie2 best-match", "ot-wg-mutated-se-mapping-best-match-sorted"],
    #                    ["Bowtie2 report-all", "ot-wg-mutated-se-mapping-report-all-sorted"],
    #                    ["MMR", "corrected-other-3mis-mmr-sorted"],
    #                    ["REMU", "corrected-ot-wg-mutated-se-mapping-remu-sorted"]]
    #
    # evaluation_results = open("./results/variants-comparison-OT-wg-70-140.txt", 'w')

    for variant_caller in variant_caller_lst:
        vcf_files = copy.deepcopy(vcf_files_names)
        for i in range(len(vcf_files)):
            vcf_files[i][1] = file_path + vcf_files[i][1] + "-variants-{}.vcf".format(variant_caller[1])

        original_variants = read_vcf_file(file_path +
                                          "mtb-h37rv-mapping-best-match-sorted-variants-{}-original.vcf".format(variant_caller[1]))
        # print(original_variants)

        comparison_output = \
            compare_variants(
                "/mnt/e/Codes/bayesian-update/data/genomes/MTB-H37Rv-back-mutated-full-mutations.txt",
                vcf_files, original_variants)

        # comparison_output = \
        #     compare_variants("/mnt/e/Codes/bayesian-update/data/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full-mutations.txt",
        #                      vcf_files)
        # comparison_output = \
        #     compare_variants("/mnt/e/Codes/bayesian-update/data/genomes/ot-whole-genome-mutated-70-140-mutations.txt",
        #                      vcf_files)
        evaluation_results.write("{}\n".format(variant_caller[0]))
        evaluation_results.write(comparison_output)

    evaluation_results.close()

# find_unique_reads("./read-mapping/mtb-whole-genome-mutated-70-140/mtb-wg-mutated-se-mapping-report-all.sam")

# find_unique_reads("./read-mapping/mtb-whole-genome-mutated-70-140/mtb-wg-mutated-se-mapping-report-all.sam")