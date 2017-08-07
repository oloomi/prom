import math

base_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def initial_counts(genome_seq, coverage):
    """
    If we see a G in reference genome, we set its initial count (prior) to [1,1,10,1]
    :param ref_genome_file:
    :return: A list of lists
    """
    base_counts = []
    for base in genome_seq:
        # The fifth one is the actual count of the base in the reference genome
        initial_base_count = [1, 1, 1, 1, coverage]
        initial_base_count[base_index[base]] = 255
        # initial_base_count = [100, 100, 100, 100]
        # initial_base_count[base_index[base]] = 700
        base_counts.append(initial_base_count)

    return base_counts


def update_counts(base_counts, selected_mapping, coverage):
    """
    Updates base counts for each position at reference
    """
    mapping_start_pos = selected_mapping[1]
    read_seq = selected_mapping[2]
    for index, base in enumerate(read_seq):
        # We treat N's as a match to the reference genome
        if base != 'N' and base_counts[mapping_start_pos + index][base_index[base]] != 255:
            if base_counts[mapping_start_pos + index][base_index[base]] < coverage:
                base_counts[mapping_start_pos + index][base_index[base]] += 1
            if base_counts[mapping_start_pos + index][4] > 1:
                base_counts[mapping_start_pos + index][4] -= 1
        else:   # If it's a match to reference
            if base_counts[mapping_start_pos + index][4] < coverage:
                base_counts[mapping_start_pos + index][4] += 1


    return True


def calc_log_mapping_prob(base_counts, mapping_start_pos, read_seq):
    """
    Calculates mapping probability for one read to one location
    It returns the sum of log of probabilities for each base
    """
    log_mapping_prob = 0
    for index, base in enumerate(read_seq):
        # If there is a match with the reference genome at this position
        if base == 'N' or base_counts[mapping_start_pos + index][base_index[base]] == 255:
            base_prob = base_counts[mapping_start_pos + index][4] / 256     # 1000
        else:
            base_prob = base_counts[mapping_start_pos + index][base_index[base]] / 256
        # base_prob = base_counts[mapping_start_pos + index][base_index[base]] / \
        #             sum(base_counts[mapping_start_pos + index])
        log_mapping_prob += math.log(base_prob)

    return log_mapping_prob
