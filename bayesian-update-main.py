from collections import defaultdict
import random
from math import log, exp

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
        initial_base_count = [1, 1, 1, 1]
        initial_base_count[base_index[base]] = 10
        base_counts.append(initial_base_count)

    return base_counts


def read_sam_file(sam_file_name):
    # A dictionary like: {read_id: [list of mappings]}
    read_alignments_dict = defaultdict(list)

    total_reads = 0
    not_mapped_reads = 0  # Number of reads that are not mapped to any location
    invalid_cigars = 0  # Number of read mapping locations that contain indels, etc.

    # Reading the SAM file and creating a dictionary of read_id : alignment
    with open(sam_file_name) as sam_file:
        for line in sam_file:
            # Skip header lines
            if line[0] == "@":
                continue

            total_reads += 1
            fields = line.rstrip().split("\t")
            read_id = fields[0]  # QNAME: Query template NAME
            cigar = fields[5]  # CIGAR string (ie. alignment)
            pos = int(fields[3])  # 1-based leftmost mapping POSition
            md_z = fields[-2]  # Alignment
            read_seq = fields[9]    # Read sequence
            # * means no alignment for a read
            if cigar != "*":
                if cigar == '150M':
                    # Store all alignments of a read
                    read_alignments_dict[read_id].append((pos, cigar, md_z, read_seq))
                else:
                    # print("Invalid CIGAR:", cigar)
                    invalid_cigars += 1
            else:
                not_mapped_reads += 1
    print("Total number of reads:", total_reads)
    print("Number of reads with non 150M CIGAR:", invalid_cigars)
    print("Number of reads not mapped:", not_mapped_reads)
    print("Number of reads in use:", len(read_alignments_dict))
    return read_alignments_dict


def calc_log_mapping_prob(base_counts, mapping_start_pos, read_seq):
    """
    Calculates mapping probability for one read to one location
    It returns the sum of log of probabilities for each base
    """
    log_mapping_prob = 0
    for index, base in enumerate(read_seq):
        base_prob = base_counts[mapping_start_pos + index][base_index[base]] / \
                    sum(base_counts[mapping_start_pos + index])
        log_mapping_prob += log(base_prob)

    return log_mapping_prob

def select_mapping():
    return True

def bayesian_update(ref_genome_file, sam_file):
    base_counts = initial_counts(ref_genome_file)
    reads_dict = read_sam_file(sam_file)
    multi_read_probs = defaultdict(list)

    # Updating the prior (initial counts) for uniquely mapped reads
    unique_reads = []
    for read_id, mappings in reads_dict.items():
        if len(mappings) == 1:
            unique_reads.append(read_id)
            mapping_start_pos = mappings[0][0] - 1
            read_seq = mappings[0][3]
            for pos_in_read, base in enumerate(read_seq):
                # We find the base counts for that position in reference genome and
                # we update it by incrementing the count for the base in the read
                base_counts[mapping_start_pos + pos_in_read][base_index[base]] += 1

    # Removing uniquely mapped reads
    for read_id in unique_reads:
        del reads_dict[read_id]

    multi_reads = sorted(reads_dict.keys())
    print("Number of multi-reads:", len(multi_reads))
    random.seed(123)

    # For each multi-read selected by random
    for i in range(10):
        read_id = random.choice(multi_reads)
        # For each of its mapping location, we calculate the posterior probability
        mapping_probs = []
        for mapping in reads_dict[read_id]:
            (mapping_start_pos, read_seq) = (mapping[0] - 1, mapping[3])
            mapping_probs.append((calc_log_mapping_prob(base_counts, mapping_start_pos, read_seq), mapping_start_pos))

        print(mapping_probs)

    print("Done!")
    return True


bayesian_update("./read-mapping/mtb-genome-extract.fna", "./read-mapping/mtb-single-end-mapping-report-all.sam")

