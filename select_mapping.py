from math_stats import *
from operator import itemgetter
import random


def select_mapping(mappings):
    """
    Normalises mapping log probabilities and selects one mapping location stochastically
    Note: it modifies the logs in the main list and converts them to normalised probabilities
    """
    # Sort mappings based on likelihood
    mappings.sort(key=itemgetter(-1))
    # Calculating log sum
    log_probs = [mapping[-1] for mapping in mappings]
    log_probs_total = logSum(log_probs)
    for mapping in mappings:
        # Normalizing log likelihood
        mapping[-1] = math.exp(mapping[-1] - log_probs_total)

    rand_num = random.uniform(0, 1)
    for mapping in mappings[:-1]:
        if rand_num < mapping[-1]:
            return mapping
    return mappings[-1]


def select_final_mapping(mappings, prob_threshold):
    """
    Returns the most probable mapping location after all runs
    """
    # Sort mappings based on likelihood
    mappings.sort(key=itemgetter(-2), reverse=True)
    selected_mapping = mappings[0]
    last_tie_index = 0

    # All max probabilities appear at the beginning of the list; find the index of last one
    for mapping in mappings[1:]:
        # New idea!
        # If it has the same probability as highest probability
        if mapping[-2] > (selected_mapping[-2] - prob_threshold):
            last_tie_index += 1
        else:
            break
    # If there are multiple locations with the same probability
    if last_tie_index > 0:
        selected_mapping = mappings[random.randrange(last_tie_index + 1)]
    return selected_mapping
