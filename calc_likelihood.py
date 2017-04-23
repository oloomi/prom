import math


base_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def initial_counts(ref_genome_file):
    """
    If we see a G in reference genome, we set its initial count (prior) to [1,1,10,1]
    :param ref_genome_file:
    :return: A list of lists
    """
    genome_seq = ""
    with open(ref_genome_file) as ref_genome:
        for line in ref_genome:
            # Skip header lines
            if line[0] == ">":
                continue
            else:
                genome_seq += line.rstrip()

    base_counts = []
    for base in genome_seq:
        # initial_base_count = [1, 1, 1, 1]
        # initial_base_count[base_index[base]] = 2
        initial_base_count = [100, 100, 100, 100]
        # initial_base_count[base_index[base]] = 10
        initial_base_count[base_index[base]] = 700
        base_counts.append(initial_base_count)

    return base_counts


def update_counts(base_counts, selected_mapping):
    """
    Updates base counts for each position at reference
    """
    mapping_start_pos = selected_mapping[1]
    read_seq = selected_mapping[2]
    for index, base in enumerate(read_seq):
        # base_counts[mapping_start_pos + index][base_index[base]] += 5
        base_counts[mapping_start_pos + index][base_index[base]] += 1
    return True


def calc_log_mapping_prob(base_counts, mapping_start_pos, read_seq):
    """
    Calculates mapping probability for one read to one location
    It returns the sum of log of probabilities for each base
    """
    log_mapping_prob = 0
    for index, base in enumerate(read_seq):
        base_prob = base_counts[mapping_start_pos + index][base_index[base]] / \
                    sum(base_counts[mapping_start_pos + index])
        log_mapping_prob += math.log(base_prob)

    return log_mapping_prob
