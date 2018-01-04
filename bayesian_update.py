import copy
from array import array
from calc_likelihood import *
from sam_file import *
from select_mapping import *


def bayesian_update(ref_genome_file, sam_file, output_file, num_runs, prob_threshold):
    """
    Resolve multimappings using Bayesian updating
    :param ref_genome_file: the path to the FASTA file of the reference genome
    :param sam_file: the input SAM file containing all candidate mappings for multireads and unique reads
    :param output_file: the path to the output SAM file with correct mappings
    :param prob_threshold: the acceptable difference in probabilities of equally good mapping locations
    :param num_runs: number of runs
    :return: True on normal exit
    """
    # 0. Reading reference genome FASTA file and creating pseudo-count vectors according to the size of genome
    print("Reading reference genome...")
    genome_header, genome_seq = read_genome(ref_genome_file)
    initial_base_counts = create_count_arrays(len(genome_seq))

    # 1. Reading the SAM file containing all mappings and
    # 1.1 Writing unique mappings to output file and updating pseudo-counts accordingly
    # 1.2 Identifying multireads and creating the dictionary of multimappings
    print("Reading all mappings...")
    multireads_dict, read_len, read_counts = read_all_mappings(sam_file, genome_seq, output_file,
                                                               initial_base_counts)

    print("Number of reads with non {}M CIGAR: {}".format(read_len, len(read_counts['unsupported'])))
    print("Number of reads not mapped:", len(read_counts['unmapped']))
    total_reads = read_counts['unique'] + len(multireads_dict)
    print("Number of reads in use:", total_reads)
    print(" {} unique reads".format(read_counts['unique']))
    print(" {} multireads with {} candidate mappings".format(len(multireads_dict), read_counts['multi']))

    # Estimated average depth of coverage
    coverage = int(total_reads * read_len / len(genome_seq))
    print("Average depth of coverage according to mapped reads: {}".format(coverage))

    # 1.3 Correct initial counts based on coverage
    process_initial_counts(initial_base_counts, coverage, genome_seq)

    # 1.4 Check for edit distance of candidate mappings of each multiread
    initially_resolved_multireads = []
    with open(output_file, 'a') as out_file:
        for read_id, mappings in multireads_dict.items():
            mappings_filtered = filter_alignments(mappings, 3)
            multireads_dict[read_id] = mappings_filtered
            # By removing low quality alignments, it has turned to a unique mapping
            if len(mappings_filtered) == 1:
                # We treat it like a unique read and update the counts
                update_counts(initial_base_counts, mappings_filtered[0], coverage, genome_seq)
                sam_fields = mappings_filtered[0][:-2]
                out_file.write("\t".join(sam_fields[0:sam_col['pos']] + [str(sam_fields[sam_col['pos']] + 1)] +
                                         sam_fields[sam_col['pos'] + 1:]))
                initially_resolved_multireads.append(read_id)

    for read_id in initially_resolved_multireads:
        del multireads_dict[read_id]

    print("Number of initially resolved multireads: {}".format(len(initially_resolved_multireads)))
    # return True
    random.seed(12)
    multi_reads = sorted(multireads_dict.keys())

    # 2. Multimapping resolution
    print("Resolving multimappings...")
    for run_number in range(num_runs):
        print("Run # {} ...".format(run_number + 1))
        random.seed(run_number)
        # Reinitialising the pseudo-counts
        base_counts = copy.deepcopy(initial_base_counts)
        # Shuffling the list of multireads
        random.shuffle(multi_reads)

        # Iterations
        num_iterations = len(multi_reads) * 1
        for i in range(num_iterations):
            # For a multi-read selected by random from the set of multi-reads
            # read_id = random.choice(multi_reads)
            read_id = multi_reads[i]

            # For each of its mapping location, we calculate the mapping probability
            for mapping in multireads_dict[read_id]:
                # Probability is saved in mapping[-1]
                calc_log_mapping_prob(base_counts, mapping, coverage, genome_seq)

            # Selecting one location stochastically
            # After select_mapping, each mapping[-1] is modified and contains probabilities instead of log-probs
            selected_mapping = select_mapping(multireads_dict[read_id])

            # Updating base counts for selected location
            update_counts(base_counts, selected_mapping, coverage, genome_seq)

        # After all iterations are done
        # Recalculate the probabilities for each multi-read and each of it's mapping locations based on latest counts
        for read_id in multi_reads:
            # For each of its mapping locations, we find the mapping probability
            for mapping in multireads_dict[read_id]:
                # Probability is saved in mapping[-1]
                calc_log_mapping_prob(base_counts, mapping, coverage, genome_seq)

            # Just run select_mapping to normalize probabilities
            select_mapping(multireads_dict[read_id])

            # Add normalised probabilities to sum of probabilities in all iterations
            for mapping in multireads_dict[read_id]:
                mapping[-2] += mapping[-1]


    # 3. Finding average probability over runs and selecting one location
    # After all runs are done
    # For each multi-read, find the average probabilities for each mapping location
    with open(output_file, 'a') as out_file:
        for read_id, mappings in multireads_dict.items():
            sum_probs = 0
            for mapping in mappings:
                # average probability over all runs for this location
                mapping[-2] /= num_runs
                sum_probs += mapping[-2]
            # Normalising average probabilities over all candidate locations for this multiread
            for mapping in mappings:
                mapping[-2] /= sum_probs

            if read_id == "gi|448814763|ref|NC_000962.3|MTB|-1622":
                probs = [(m[sam_col['pos']], m[-2]) for m in mappings]
                print(sorted(probs))
                # counts_print(initial_base_counts, 2110, 2160)
                # print('')
                # counts_print(base_counts, 2110, 2160)

            # Selecting final mapping
            final_mapping = select_final_mapping(mappings, prob_threshold)
            # Add 1 to position and write to file
            out_file.write("\t".join(final_mapping[0:sam_col['pos']] + [str(final_mapping[sam_col['pos']] + 1)] +
                                     final_mapping[sam_col['pos'] + 1:-2]))
    return True
