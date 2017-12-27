import copy
from array import array
from calc_likelihood import *
from sam_file import *
from select_mapping import *


def bayesian_update(ref_genome_file, sam_file, output_file):
    """
    Resolve multimappings using Bayesian updating
    :param ref_genome_file: the path to the FASTA file of the reference genome
    :param sam_file: the input SAM file containing all candidate mappings for multireads and unique reads
    :param output_file: the path to the output SAM file with correct mappings
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

    # 1.3 Check initial base counts
    process_initial_counts(initial_base_counts, coverage)

    return True

    # 1. Finding initial counts
    # initial_base_counts = initial_counts(genome_seq, coverage)

    multi_reads_final_location = defaultdict(int)

    # Updating the prior (initial counts) for uniquely mapped reads
    # New: and also filtering may make a multi-read a unique read
    unique_reads = []
    initially_resolved_multireads = []
    random.seed(12)

    for read_id, mappings in reads_dict.items():
        if len(mappings) == 1:
            unique_reads.append(read_id)
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
                initially_resolved_multireads.append(read_id)
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

    # {read_id: {run_number: [probabilities]}, ...}
    # multi_read_probs = defaultdict(lambda: defaultdict(list))
    # [[(pos, prob),...], ... ] for run 1, 2, 3, etc.
    multi_read_probs = defaultdict(list)
    # random_seeds = [12, "Hi", 110, "Bye", 1, 33, 5, 14, 313, 777]

    # Logging counts in each run
    counts_log_file = open("{}-log-counts.txt".format(output_file[:-4]), 'w')

    # Trying new idea
    sum_pr = 0
    cnt_pr = 0

    # 2. Sampling
    # 10 Runs with 5000 iterations in each
    for run_number in range(10):
        print("Run # {} ...".format(run_number))
        # random.seed(random_seeds[run_number])
        random.seed(run_number)
        base_counts = copy.deepcopy(initial_base_counts)

        # Trying new idea
        random.shuffle(multi_reads)

        # Iterations
        num_iterations = len(multi_reads) * 1
        for i in range(num_iterations):
            # For a multi-read selected by random from the set of multi-reads
            # read_id = random.choice(multi_reads)
            read_id = multi_reads[i]

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
            selected_mapping = select_mapping(mapping_probs)

            # test
            # if read_id == 'gi|448814763|ref|NC_000962.3|MTB|-673':
            #     sum_pr += mapping_probs[1][0]
            #     cnt_pr += 1
            #     print(mapping_probs[1][0], sum_pr / cnt_pr)

            # Updating base counts for selected location
            update_counts(base_counts, selected_mapping, coverage)

        # After all iterations are done
        # Save the probabilities for each multi-read and each of it's mapping locations
        for read_id in multi_reads:
            # For each of its mapping location, we find the mapping probability
            mapping_probs = []  # [(position, probability), ...]
            for mapping in reads_dict[read_id]:
                (mapping_start_pos, read_seq) = (mapping[0] - 1, mapping[3])
                mapping_probs.append([calc_log_mapping_prob(base_counts, mapping_start_pos, read_seq),
                                      mapping_start_pos, read_seq])
            # Just run select_mapping to normalize probabilities
            select_mapping(mapping_probs)
            pos_prob_list = [(mp[1], mp[0]) for mp in mapping_probs]
            # Sorted by position
            multi_read_probs[read_id].append(sorted(pos_prob_list))

        # Log: writing counts to file for logging
        counts_log_file.write("Run {}\n".format(run_number + 1))
        for pos, base_vector in enumerate(base_counts):
            counts_log_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(pos+1, base_vector[0], base_vector[1],
                                                                    base_vector[2], base_vector[3], base_vector[4]))


    # 3. Finding average probability over runs and selecting one location
    # After all runs are done
    final_multiread_probs = defaultdict(list)
    # For each multi-read, find the average probabilities for each mapping location
    for read_id, mapping_probs_list in multi_read_probs.items():
        # mapping_probs_list
        # [[run_1], [run_2], [run_3], ...]
        # [[(pos, prob), (pos, prob), ...], [run_2], ...]
        avg_probs = [0 for _ in range(len(mapping_probs_list[0]))]
        # Sum of probabilities at each position over runs
        for mapping_probs in mapping_probs_list:
            for i, pos_prob in enumerate(mapping_probs):
                avg_probs[i] += pos_prob[1]
        # Average at each position over all runs
        for i in range(len(avg_probs)):
            avg_probs[i] /= len(mapping_probs_list)
        # Normalizes probabilities found after averaging over runs
        sum_probs = sum(avg_probs)
        for i, pos_prob in enumerate(mapping_probs_list[0]):
            # (normalized probability, position)
            avg_probs[i] = (avg_probs[i] / sum_probs, pos_prob[0])
        final_multiread_probs[read_id] = avg_probs

    # The final selected mapping location for reach multi-read and write the log in file
    with open("{}-log.txt".format(output_file[:-4]), 'w') as log_file:
        # multi_reads_final_location = defaultdict(int)
        for read_id, mapping_probs in final_multiread_probs.items():
            if read_id == 'gi|448814763|ref|NC_000962.3|MTB|-1622':
                probs = [prob[0] for prob in mapping_probs]
                print(mapping_probs)
                print(probs)
            # Selecting final mapping location
            best_mapping_location = select_final_mapping(mapping_probs)
            if read_id == 'gi|448814763|ref|NC_000962.3|MTB|-1622':
                print(best_mapping_location)
            # best_mapping_location = select_final_mapping_stochastic(mapping_probs)
            multi_reads_final_location[read_id] = best_mapping_location[1] + 1

            # Writing log to file
            # loc_prob = [(loc, round(prob, 2)) for prob, loc in final_multiread_probs[read_id]]
            # log_file.write("{}\t{}\t{}\t{}\n".format(read_id, multi_reads_final_location[read_id], 'M', loc_prob))
            # for i, run_log in enumerate(multi_read_probs[read_id]):
            #     loc_prob_logs = [(loc, round(prob, 2)) for loc, prob in run_log]
            #     log_file.write("{}\t{}\t{}\t{}\n".format(read_id, multi_reads_final_location[read_id], i,
            #                                              loc_prob_logs))

    # Writing final results to a SAM file
    write_sam_file(multi_reads_final_location, sam_file, output_file)

    return True
