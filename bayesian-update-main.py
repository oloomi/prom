from collections import defaultdict
from samfile import *
import random
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
        initial_base_count = [1, 1, 1, 1]
        initial_base_count[base_index[base]] = 10
        base_counts.append(initial_base_count)

    return base_counts


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


def logAdd(a, b):
    """
    return log(exp(a) + exp(b))
    Source: https://github.com/drtconway/pykmer
    """
    x = max(a, b)
    y = min(a, b)
    w = y - x
    return x+math.log1p(math.exp(w))


def logSum(xs):
    """
    return log(sum([exp(x) for x in xs]))
    Source: https://github.com/drtconway/pykmer
    """
    assert len(xs) > 0
    y = xs[0]
    for x in xs[1:]:
        y = logAdd(y, x)
    return y


def select_mapping(mapping_probs):
    """
    Normalises mapping log probabilities and selects one mapping location stochasticaly
    Note: it modifies the logs in the main list
    """
    mapping_probs.sort()
    log_probs = [mapping[0] for mapping in mapping_probs]
    log_probs_total = logSum(log_probs)
    for mapping in mapping_probs:
        # Normalizing log likelihood
        mapping[0] = math.exp(mapping[0] - log_probs_total)

    rand_num = random.uniform(0, 1)
    for mapping in mapping_probs[:-1]:
        if rand_num < mapping[0]:
            return mapping
    return mapping_probs[-1]


def best_mapping(mapping_probs):
    """
     Selecting the best mapping among candidate locations for a multi-read after last iteration
    """
    mapping_probs = sorted(mapping_probs, reverse=True)
    selected_prob = mapping_probs[0]
    last_tie_index = 0

    # All max probabilities appear at the beginning of the list; find the index of last one
    for mapping in mapping_probs[1:]:
        # If it has the same probability as highest probability
        if mapping[0] == selected_prob[0]:
            last_tie_index += 1
        else:
            break
    # If there are multiple locations with the same probability
    if last_tie_index > 0:
        selected_prob = mapping_probs[random.randrange(last_tie_index + 1)]

    # For probability normalisation
    log_probs = [mapping[0] for mapping in mapping_probs]
    log_probs_total = logSum(log_probs)

    selected_prob[0] = math.exp(selected_prob[0] - log_probs_total)

    return selected_prob


def update_counts(base_counts, selected_mapping):
    """
    Updates base counts for each position at reference
    """
    mapping_start_pos = selected_mapping[1]
    read_seq = selected_mapping[2]
    for index, base in enumerate(read_seq):
        base_counts[mapping_start_pos + index][base_index[base]] += 1
    return True


def bayesian_update(ref_genome_file, sam_file, output_file):
    """
    Assigning a multi-read to a mapping location using Bayesian update
    :param ref_genome_file:
    :param sam_file:
    :return:
    """
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
    for i in range(5000):
        read_id = random.choice(multi_reads)
        # For each of its mapping location, we calculate the posterior probability
        mapping_probs = []  # (probability, position, read_seq)
        for mapping in reads_dict[read_id]:
            (mapping_start_pos, read_seq) = (mapping[0] - 1, mapping[3])
            mapping_probs.append([calc_log_mapping_prob(base_counts, mapping_start_pos, read_seq),
                                  mapping_start_pos, read_seq])

        # Selecting one location statistically
        selected_mapping = select_mapping(mapping_probs)

        # For tracking convergence: latest probabilities, first
        multi_read_probs[read_id].append([mapping[0] for mapping in mapping_probs])

        # Updating base counts for selected location
        update_counts(base_counts, selected_mapping)

    with open(output_file, 'w') as out_file:
        for read_id, normalized_probs in multi_read_probs.items():
            out_file.write("*** {} ***\n\n".format(read_id))
            for prob in normalized_probs:
                out_file.write("{}\n".format(prob))
            out_file.write("\n")

    # Find the final mapping location
    multi_reads_final_location = defaultdict(int)
    for read_id in multi_reads:
        # For each of its mapping location, we find the mapping probability
        mapping_probs = []  # (probability, position, read_seq)
        for mapping in reads_dict[read_id]:
            (mapping_start_pos, read_seq) = (mapping[0] - 1, mapping[3])
            mapping_probs.append([calc_log_mapping_prob(base_counts, mapping_start_pos, read_seq),
                                  mapping_start_pos, read_seq])

        # Selecting one location statistically
        best_mapping_location = best_mapping(mapping_probs)
        multi_reads_final_location[read_id] = best_mapping_location[1] + 1

        print("_______________________")
        print(read_id, best_mapping_location[1] + 1, round(best_mapping_location[0], 4))
        print(".......................")
        for mapping in reads_dict[read_id]:
            print(mapping[0:3])

    # Writing final results to a SAM file
    write_sam_file(multi_reads_final_location, sam_file, "corrected-mappings-mtb.sam")

    return True


bayesian_update("./read-mapping/mtb-genome-extract.fna", "./read-mapping/mtb-single-end-mapping-report-all.sam",
                "mtb-normalised-probs.txt")

