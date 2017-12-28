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
    """
    Updates the initial counts based on nucleotides in the unique reads
    """
    mapping_start_pos = selected_mapping[sam_col['pos']]
    read_seq = selected_mapping[sam_col['seq']]
    for index, base in enumerate(read_seq):
        ref_base = genome_seq[mapping_start_pos + index]
        # We treat N's as a match to the reference genome
        if base != 'N':
            base_counts[base_index[base]][mapping_start_pos + index] += 1
            base_counts[base_index[ref_base]][mapping_start_pos + index] -= 1
        else:   # Update reference base
            base_counts[base_index[ref_base]][mapping_start_pos + index] += 1
    return True


def process_initial_counts(base_counts, coverage):
    """
    Checks the initial counts and corrects counts greater than coverage or less than one
    """
    genome_size = len(base_counts[0])
    for i in range(4):
        for j in range(genome_size):
            if base_counts[i][j] > coverage:
                base_counts[i][j] = coverage
            elif base_counts[i][j] < 1:
                base_counts[i][j] = 1
    return True


def update_counts(base_counts, selected_mapping, coverage, genome_seq):
    """
    Updates base counts for each position at reference according to the alignment of multiread
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


def calc_log_mapping_prob(base_counts, mapping, coverage, genome_seq):
    """
    Calculates mapping probability for one read to one location
    It returns the sum of log of probabilities for each base
    """
    log_mapping_prob = 0
    for index, base in enumerate(mapping[sam_col['seq']]):
        if base != 'N':
            base_prob = base_counts[base_index[base]][mapping[sam_col['pos']] + index] / (coverage + 1)
        else:
            ref_base = genome_seq[mapping[sam_col['pos']] + index]
            base_prob = base_counts[base_index[ref_base]][mapping[sam_col['pos']] + index] / (coverage + 1)

        log_mapping_prob += math.log(base_prob)

    mapping[-1] = log_mapping_prob

    return True
