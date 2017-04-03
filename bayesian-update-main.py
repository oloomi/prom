from collections import defaultdict
from samfile import *
import random
import math
import copy

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
        base_counts[mapping_start_pos + index][base_index[base]] += 5
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


def select_final_mapping(mapping_probs):
    """
    Returns the most probable mapping location after all runs
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
    return selected_prob


def bayesian_update(ref_genome_file, sam_file, output_file):
    """
    Assigning a multi-read to a mapping location using Bayesian update
    :param ref_genome_file:
    :param sam_file:
    :return:
    """
    initial_base_counts = initial_counts(ref_genome_file)
    reads_dict = read_sam_file(sam_file)

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
                initial_base_counts[mapping_start_pos + pos_in_read][base_index[base]] += 1

    # Removing uniquely mapped reads
    for read_id in unique_reads:
        del reads_dict[read_id]

    multi_reads = sorted(reads_dict.keys())
    print("Number of multi-reads:", len(multi_reads))


    # {read_id: {run_number: [probabilities]}, ...}
    # multi_read_probs = defaultdict(lambda: defaultdict(list))
    # [[(pos, prob),...], ... ] for run 1, 2, 3, etc.
    multi_read_probs = defaultdict(list)
    random_seeds = [12, "Hi", 110, "Bye", 1, 33, 5, 14, 313, 777]

    # 3 Runs with 5000 iterations in each
    for run_number in range(10):
        # random.seed(random_seeds[run_number])
        random.seed(run_number)
        base_counts = copy.deepcopy(initial_base_counts)

        # Iterations
        for i in range(5000):
            # For a multi-read selected by random from the set of multi-reads
            read_id = random.choice(multi_reads)

            # !!! For debugging !!!
            # if read_id == 'gi|20000|ref|NC_000962.3|3930000_3949999|-497':
            #     print('\n', run_number, i)

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

            # !!! For debugging !!!
            # if read_id == 'gi|20000|ref|NC_000962.3|3930000_3949999|-497':
            #     print([(mp[0], mp[1])for mp in mapping_probs])
            #     print(base_counts[3641], base_counts[4244])

            # For tracking convergence: latest probabilities, first
            # multi_read_probs[read_id][run_number].append([mapping[0] for mapping in mapping_probs])

            # Updating base counts for selected location
            update_counts(base_counts, selected_mapping)

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

    # After all runs are done
    final_multiread_probs = defaultdict(list)
    # For each multi-read, find the average probabilities for each mapping location
    for read_id, mapping_probs_list in multi_read_probs.items():
        # mapping_probs_list
        # [[run_1], [run_2], [run_3], ...]
        # [[(pos, prob), (pos, prob), ...], [run_2], ...]
        avg_probs = [0 for i in range(len(mapping_probs_list[0]))]
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

    # The final selected mapping location for reach multi-read
    multi_reads_final_location = defaultdict(int)

    for read_id, mapping_probs in final_multiread_probs.items():
        best_mapping_location = select_final_mapping(mapping_probs)
        multi_reads_final_location[read_id] = best_mapping_location[1] + 1

        print("________________")
        print(read_id)
        print(multi_reads_final_location[read_id])
        print("................")
        print(final_multiread_probs[read_id])
        print("................")
        for run_log in multi_read_probs[read_id]:
            print(run_log)

    # Writing final results to a SAM file
    write_sam_file(multi_reads_final_location, sam_file, output_file)

    return True


# bayesian_update("./data/genomes/mtb-genome-extract.fna",
#                 "./read-mapping/mtb-normal/mtb-normal-se-mapping-report-all.sam",
#                 "./read-mapping/mtb-normal/corrected-mappings-mtb-normal-700-100-5.sam")

bayesian_update("./data/genomes/mtb-genome-extract-mutated.fna",
                "./read-mapping/mtb-mutated/mtb-mutated-se-mapping-report-all.sam",
                "./read-mapping/mtb-mutated/corrected-mappings-mtb-mutated-700-100-5.sam")

# unique_reads("./read-mapping/mtb-mutated/mtb-mutated-se-mapping-report-all.sam")
