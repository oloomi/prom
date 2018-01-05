import math
from array import array
from copy import deepcopy
from settings import *


def create_count_arrays(genome_seq):
    """
    Returns a dictionary of lists that contains four arrays corresponding to the pseudo-counts of the four nucleotides
    along each chromosome in the genome
    Format: {chromosome_name: [[A counts], [C counts], [G counts], [T counts]} ie. chr:nucleotide:pos
    """
    base_counts = {}
    for name, seq in genome_seq.items():
        one_array = array('i', [0] * len(seq))
        base_counts[name] = [deepcopy(one_array), deepcopy(one_array), deepcopy(one_array), deepcopy(one_array)]
    return base_counts


def initial_counts(base_counts, selected_mapping, genome_seq):
    """
    Updates the initial counts based on nucleotides in the unique reads
    """
    chrom = selected_mapping[sam_col['rname']]
    mapping_start_pos = selected_mapping[sam_col['pos']]
    read_seq = selected_mapping[sam_col['seq']]
    base_qual = selected_mapping[sam_col['qual']]
    for index, base in enumerate(read_seq):
        ref_base = genome_seq[chrom][mapping_start_pos + index]
        # We treat N's and low quality score bases as a match to the reference genome
        # Phred-scale quality score
        if (ord(base_qual[index]) - 33) < base_qual_threshold or base == 'N':
            base_counts[chrom][base_index[ref_base]][mapping_start_pos + index] += 1
        else:
            base_counts[chrom][base_index[base]][mapping_start_pos + index] += 1
            # base_counts[base_index[ref_base]][mapping_start_pos + index] -= 1
    return True


def process_initial_counts(base_counts, coverage, genome_seq):
    """
    Checks the initial counts and corrects counts greater than coverage or less than one
    """
    for chrom in genome_seq.keys():
        for i in range(len(genome_seq[chrom])):
            nucleotides = ['A', 'C', 'G', 'T']
            ref_base = genome_seq[chrom][i]
            nucleotides.remove(ref_base)
            sum_alt_counts = 0
            # Sum of counts for alternative bases
            for base in nucleotides:
                sum_alt_counts += base_counts[chrom][base_index[base]][i]
            ref_count = coverage + base_counts[chrom][base_index[ref_base]][i] - sum_alt_counts
            base_counts[chrom][base_index[ref_base]][i] = ref_count
            for b in range(4):
                if base_counts[chrom][b][i] > coverage:
                    base_counts[chrom][b][i] = coverage
                elif base_counts[chrom][b][i] < 1:
                    base_counts[chrom][b][i] = 1

    return True


def update_counts(base_counts, selected_mapping, coverage, genome_seq):
    """
    Updates base counts for each position at reference according to the alignment of multiread
    """
    chrom = selected_mapping[sam_col['rname']]
    mapping_start_pos = selected_mapping[sam_col['pos']]
    read_seq = selected_mapping[sam_col['seq']]
    base_qual = selected_mapping[sam_col['qual']]
    # mapping_start_pos = selected_mapping[1]
    # read_seq = selected_mapping[2]
    for index, base in enumerate(read_seq):
        ref_base = genome_seq[chrom][mapping_start_pos + index]
        # We treat N's and low quality score bases as a match to the reference genome
        if base == ref_base or (ord(base_qual[index]) - 33) < base_qual_threshold or base == 'N':
            if base_counts[chrom][base_index[ref_base]][mapping_start_pos + index] < coverage:
                base_counts[chrom][base_index[ref_base]][mapping_start_pos + index] += 1
        else:
            if base_counts[chrom][base_index[base]][mapping_start_pos + index] < coverage:
                base_counts[chrom][base_index[base]][mapping_start_pos + index] += 1
            # Update the count for the reference base
            if base_counts[chrom][base_index[ref_base]][mapping_start_pos + index] > 1:
                base_counts[chrom][base_index[ref_base]][mapping_start_pos + index] -= 1
    return True


def counts_print(base_counts, chrom, start_pos, end_pos):
    offsets = [str(i) for i in range(end_pos - start_pos + 1)]
    print('\t'.join(offsets))
    for i in range(4):
        # print(base_counts[i][start_pos - 1:  end_pos])
        print('\t'.join([str(c) for c in base_counts[chrom][i][start_pos - 1:  end_pos]]))


def calc_log_mapping_prob(base_counts, mapping, coverage, genome_seq):
    """
    Calculates mapping probability for one read to one location
    It returns the sum of log probabilities for bases along the alignment
    """
    log_mapping_prob = 0
    chrom = mapping[sam_col['rname']]
    base_qual = mapping[sam_col['qual']]
    for index, base in enumerate(mapping[sam_col['seq']]):
        if (ord(base_qual[index]) - 33) < base_qual_threshold or base == 'N':
            ref_base = genome_seq[chrom][mapping[sam_col['pos']] + index]
            base_prob = base_counts[chrom][base_index[ref_base]][mapping[sam_col['pos']] + index] / (coverage + 1)
        else:
            base_prob = base_counts[chrom][base_index[base]][mapping[sam_col['pos']] + index] / (coverage + 1)

        log_mapping_prob += math.log(base_prob)

    mapping[-1] = log_mapping_prob

    return True
