import math
from array import array
from copy import deepcopy

sam_col = {'qname': 0, 'pos': 3, 'cigar': 5, 'seq': 9, 'qual': 10}
base_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def create_count_arrays(genome_size):
    """
    Returns a list that contains four arrays corresponding to the pseudo-counts of the four nucleotides along the genome
    """
    one_array = array('i', [1] * genome_size)
    base_counts = [deepcopy(one_array), deepcopy(one_array), deepcopy(one_array), deepcopy(one_array)]
    return base_counts


def initial_counts(base_counts, selected_mapping, genome_seq):
    mapping_start_pos = selected_mapping[sam_col['pos']]
    read_seq = selected_mapping[sam_col['seq']]
    for index, base in enumerate(read_seq):
        ref_base = genome_seq[mapping_start_pos + index]
        # We treat N's as a match to the reference genome
        if base != 'N':
            base_counts[base_index[base]][mapping_start_pos + index] += 1
        else:   # Update reference base
            base_counts[base_index[ref_base]][mapping_start_pos + index] += 1
    return True


def update_counts(base_counts, selected_mapping, coverage, genome_seq):
    """
    Updates base counts for each position at reference
    """
    mapping_start_pos = selected_mapping[sam_col['pos']]
    read_seq = selected_mapping[sam_col['seq']]
    # mapping_start_pos = selected_mapping[1]
    # read_seq = selected_mapping[2]
    for index, base in enumerate(read_seq):
        ref_base = genome_seq[mapping_start_pos + index]
        # We treat N's as a match to the reference genome
        if base != 'N' and base_counts[base_index[base]][mapping_start_pos + index] != ref_base:
            if base_counts[base_index[base]][mapping_start_pos + index] < coverage:
                base_counts[base_index[base]][mapping_start_pos + index] += 1
            # Update the count for the reference base
            if base_counts[base_index[ref_base]][mapping_start_pos + index] > 1:
                base_counts[base_index[ref_base]][mapping_start_pos + index] -= 1
        else:   # If it's a match to reference
            if base_counts[base_index[ref_base]][mapping_start_pos + index] < coverage:
                base_counts[base_index[ref_base]][mapping_start_pos + index] += 1
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
